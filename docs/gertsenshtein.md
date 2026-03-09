# Gertsenshtein Effect: Physics Background and Validation Targets

## Overview

The **Gertsenshtein effect** (1962) is the resonant conversion of gravitational waves (gravitons) into electromagnetic waves (photons) — and vice versa — in an external magnetic field. The coupling arises from the standard Einstein-Maxwell system: the spacetime metric appears in the Maxwell action, so linearizing around a background magnetic field generates cross-terms between the metric perturbation h_μν and the EM perturbation a_μ.

This effect is the simplest instance of **graviton-matter wave mixing** and serves as the foundation for extensions to torsion waves, axion-photon-graviton mixing, and Chern-Simons gravity.

## Physics

### The Einstein-Maxwell System

The full Lagrangian:

```
L = (1/κ²) R - (1/4) F_{μν} F^{μν}
```

where κ² = 16πG, R is the Ricci scalar, and F_{μν} = ∂_μ A_ν - ∂_ν A_μ is the electromagnetic field strength tensor.

### Background Configuration

- **Metric**: flat Minkowski η_{μν} = diag(-1, +1, +1, +1)
- **EM field**: background 4-potential Ā_μ generating a uniform magnetic field B₀
  - For B₀ along x-axis: Ā_y = -B₀ z (Coulomb gauge)
  - Background field strength: F̄_{yz} = ∂_z Ā_y = -B₀, i.e., B_x = B₀

### Perturbation Expansion

Both fields are expanded around their backgrounds:

```
g_{μν} = η_{μν} + ε h_{μν}     (metric perturbation — graviton)
A_μ    = Ā_μ    + ε a_μ         (EM perturbation — photon)
```

The second-order action δ²(√|g| L) contains:
- **h × h terms**: graviton self-interaction (Fierz-Pauli)
- **a × a terms**: photon kinetic energy (Maxwell)
- **h × a terms**: graviton-photon coupling via F̄ — **the Gertsenshtein effect**

The coupling coefficients are proportional to the background field strength B₀.

### Gauge Fixing

- **Metric**: Transverse-traceless (TT) gauge — reduces to physical polarization modes h_+, h_×
- **EM**: Lorenz gauge ∂_μ a^μ = 0 — eliminates unphysical degrees of freedom

### Plane-Wave Reduction

For GW propagating along z-axis with B₀ along x-axis:
- Set ∂/∂x = ∂/∂y = 0 (plane-wave ansatz)
- Reduces 3+1D tensor+vector system to effective 1+1D coupled equations
- Only h_+ polarization couples to a_y (the transverse EM mode parallel to B₀)

### Conversion Probability

**Vacuum, uniform B₀, propagation distance D = c·t:**

The coupled equations from our Lagrangian (1/κ²)R - (1/4)F² are (after TT + Lorenz
gauge, plane-wave reduction, and dividing by kinetic coefficients):

```
d²t(h_+) = ∂²z(h_+) - B₀²κ² h_+ - B₀κ² ∂z(a_y)
d²t(a_y) = ∂²z(a_y) + B₀ ∂z(h_+)
```

The coupling is asymmetric: B₀κ² on h (from dividing by 1/κ² kinetic prefactor) vs
B₀ on a (unit kinetic prefactor). The eigenmode analysis gives beat frequency
Δω = B₀κ (geometric mean of couplings), yielding:

```
P(graviton → photon) = sin²(κ B₀ D / 2)
```

This is confirmed numerically (RMS error < 0.02 across multiple B₀ values).

**Note on P&R error:** Palessandro & Rothman (2023, Eq. 26; 2024, arXiv:2405.01407)
quote P = sin²(√G · B₀ · D), which is incorrect by a factor of √(4π) = 2√π ≈ 3.54
in the argument. They fail to properly account for the 1/(16πG) kinetic prefactor
when extracting the mixing frequency from their non-canonically-normalized system.
Domcke & Garcia-Cely (2023, arXiv:2310.04150) independently criticize P&R's
linearization. Our formula is confirmed by Dandoy, Lella et al. (2024,
arXiv:2406.17853), who use canonical normalization and obtain Δ_{g,γ} = √(4πG)·B₀,
matching κB₀/2 exactly. See [gertsenshtein_formula.md](gertsenshtein_formula.md)
for the full derivation and literature comparison.

**With detuning** (plasma frequency ω_p, or effective photon mass m_γ):

```
P = f² sin²(ω_m D)

where:
  Δ = (m_γ² - m_g²) / (2ω)     (mass-squared difference / 2ω)
  μ = κ B₀ / 2                  (coupling strength)
  ω_m = √(Δ² + μ²)             (mixing frequency)
  f = μ / ω_m                   (mixing angle)
```

**Weak-field limit** (μ D ≪ 1):

```
P ≈ κ² B₀² D² / 4
```

**Oscillation length**:

```
L_osc = 2π / (κ B₀)
```

### Key Physics Regimes

| Regime | Condition | P_max | Character |
| ------ | --------- | ----- | --------- |
| Weak coupling | μ D ≪ 1 | ≈ κ²B₀²D²/4 | Quadratic growth |
| Resonant | Δ = 0 (vacuum) | sin²(μ D) | Full Rabi oscillation |
| Detuned | Δ ≫ μ | ≈ (μ/Δ)² | Suppressed, fast oscillation |
| MSW resonance | ω_p² = m_g² (at some position) | Enhanced | Level crossing |

## References

### Primary Sources

| Ref | Year | Key Contribution |
| --- | ---- | ---------------- |
| [Gertsenshtein, JETP 14:84](https://www.scirp.org/reference/referencespapers?referenceid=1912195) | 1962 | Original prediction of graviton-photon conversion in B-field |
| [Boccaletti et al., Nuovo Cimento 70B:129](https://link.springer.com/article/10.1007/BF02710177) | 1970 | Finite-region boundary conditions, integral formula |
| [Raffelt & Stodolsky, PRD 37:1237](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.37.1237) | 1988 | 2×2 mixing-matrix formalism, MSW resonance analogy |
| [Ejlli, JHEP 2020:029](https://arxiv.org/abs/2004.02714) | 2020 | First exact (non-perturbative) solution |
| [Palessandro & Rothman, PDU 40:101187](https://arxiv.org/abs/2301.02072) | 2023 | Lagrangian derivation (**contains normalization error** — see [formula comparison](gertsenshtein_formula.md)) |
| [Domcke & Garcia-Cely, arXiv:2310.04150](https://arxiv.org/abs/2310.04150) | 2023 | Criticizes P&R linearization; inhomogeneous B-field calculation |
| [**Dandoy, Lella et al., arXiv:2406.17853**](https://arxiv.org/abs/2406.17853) | **2024** | **4-component mixing, canonical graviton — confirms our formula** |
| [Palessandro, arXiv:2405.01407](https://arxiv.org/abs/2405.01407) | 2024 | Extended treatment (same normalization error as 2023) |

### Extensions and Modern Applications

| Ref | Year | Key Contribution |
| --- | ---- | ---------------- |
| [Domcke & Garcia-Cely, JCAP 05:051](https://arxiv.org/abs/2312.17636) | 2024 | Inverse Gertsenshtein as HFGW probe |
| [**Domcke, Garcia-Cely & Lee, arXiv:2507.16609**](https://arxiv.org/abs/2507.16609) | **2025** | **GW scattering on 3D B-fields; Born approx. = thin-lens Boccaletti; WKB = thick-lens** |
| [Obukhov et al., arXiv:2410.01355](https://arxiv.org/abs/2410.01355) | 2024 | Photon-torsion wave conversion in Poincaré gauge theory |
| [arXiv:2507.02362](https://arxiv.org/abs/2507.02362) | 2025 | EM-torsion coupling near black holes |

### Numerical Methods

| Ref | Year | Key Contribution |
| --- | ---- | ---------------- |
| [Berlin et al., arXiv:2405.08865](https://arxiv.org/abs/2405.08865) | 2024 | Numerical analysis of resonant axion-photon mixing — directly analogous methodology; WKB vs exact regimes |

## Validation Targets

All equations are **derived from the Lagrangian** by the TIDAL pipeline — listed here only for comparison.

### Target 1: Vacuum Conversion Probability

**Source**: Eigenmode analysis of coupled h_7/a_2 system; confirmed by Dandoy, Lella et al. (2024).

- Uniform B₀, vacuum (no plasma), massless graviton
- Plane wave at known wavenumber k propagating along z
- Measure P(t) via `tidal measure --what conversion`
- Compare: `P = sin²(κ B₀ D / 2)`
- **Acceptance**: Agreement to within 1%

### Target 2: B₀ Sweep

- Sweep B₀ from weak to strong coupling
- Plot P_final vs B₀
- Compare against analytical curve
- Check oscillation length scaling: L_osc ∝ 1/B₀

### Target 3: Detuned Conversion (Phase F1)

**Status: NOT IMPLEMENTED — requires gauge-invariant mass mechanism**

Adding a Proca-type photon mass `- omegaP2/2 * a[-mu] eta[mu, nu] a[-nu]` to the full 4D
Lagrangian generates position-dependent coefficients `B₀²κ²omegaP2 z²` on the graviton equations.
These arise from the volume element expansion `(ε/2)Tr(h) × (omegaP2/2)|Ā|²` where `|Ā|² = B₀²z²`
(since `Ā_y = -B₀z` in the plane-wave gauge). With periodic BCs, z² is discontinuous at the wrap
point — incorrect for a uniform plasma model.

**Root cause**: xPert correctly evaluates the full second-order action including the metric coupling
to the background EM energy density. In the plasma context (which is a medium effect, not a
fundamental Proca field), this coupling is unphysical. No Lagrangian reformulation within the
current 4D gauge-field framework can avoid this: any gauge representation of a transverse B-field
has `|Ā|² ∝ z²` in 1D. An effective 1+1D scalar theory (hgw, aem scalars) yields constant
coefficients but is no different from the existing `coupled_scalars` example.

**What would be needed**: A gauge-invariant effective mass that couples only to the transverse
photon polarization without coupling to the background gauge potential through the volume element.
Possible approaches: (a) a Stueckelberg mass in unitary gauge with a carefully chosen gauge field,
or (b) a medium/refractive-index formulation at the level of the equations of motion (bypassing
the Lagrangian for the plasma sector). See `.github-issues-pending.md` for the tracking issue.

- **Reference**: Raffelt & Stodolsky (1988, PRD 37:1237); Dandoy, Lella et al. (2024, arXiv:2406.17853)

### Target 4: Finite-Region Scattering (Phase F2)

**Status: IMPLEMENTED** — see `examples/gertsenshtein/theory_localized.toml`

**Setup**: Gaussian B₀(z) = Bpeak × exp(-z²/(2R²)) along x-axis. The gauge potential
`Ā_y = -Bpeak × R × √(π/2) × Erf(z/(√2·R))` satisfies ∂_z(Ā_y) = -B₀(z). The
`Erf` function is supported by the TIDAL evaluator (`_eval_utils.py` `_FUNCTION_MAP`).
Constant named `Bpeak` (not `B0_peak`) — in Wolfram, `X_Y` parses as `Pattern[X, Blank[Y]]`
(a pattern expression, not a symbol), silently corrupting the symbolic computation.

**Analytical target (Boccaletti 1970, weak-field Gaussian)**:
```
P(graviton → photon) = sin²(κ/2 × ∫B₀(z)dz) = sin²(κ × Bpeak × R × √(π/2))
```

**Simulation parameters**:
- Domain: [-100, 100], 1024 points, non-periodic (Neumann BCs) — periodic BCs are
  incorrect since Ā_y is not periodic
- IC: Gaussian wavepacket (h_7) at z = -50, σ=5, k = 2.0, t_end = 120
  (stops after one pass through B-field, before wavepacket hits right boundary at t≈150)
- Measurement: `--what conversion --source h_7 --target a_2`
  (NOT `--what energy` — position-dependent coefficients cause known limitation)
- **Bounds/center require `=` syntax**: `--bounds="-100:100"`, `--ic-center=-50.0`
  (shell interprets leading `-` as flag otherwise)

**Confirmed numerical result** (defaults: κ=1, Bpeak=0.1, R=5=σ):
```
P_numerical = 0.3436  (at t=80.4, peak of single pass)
P_Boccaletti = sin²(1.0 × 0.1 × 5.0 × √(π/2)) = sin²(0.627) = 0.3432
Agreement: 0.04% — excellent
```

**2D sweep (Bpeak × R, 48 points)** — Boccaletti formula confirmed across all regimes (σ=5):

| R/σ | R | Mean err vs Boccaletti | Max err |
|-----|---|----------------------|---------|
| 0.40 | 2.0 | 0.00004 | 0.00008 |
| 0.92 | 4.6 | 0.00029 | 0.00045 |
| 1.44 | 7.2 | 0.00078 | 0.00172 |
| 1.96 | 9.8 | 0.00113 | 0.00211 |
| 2.48 | 12.4 | 0.00091 | 0.00144 |
| 3.00 | 15.0 | 0.00058 | 0.00160 |

The Boccaletti formula P = sin²(κ/2 × ∫B₀ dz) is **exact for massless vacuum conversion**
for any R/σ ratio. For massless graviton-photon conversion (no plasma), both modes travel at
c with identical dispersion (k_h = k_γ), so Δk = 0 exactly. The conversion amplitude at
q = Δk = 0 is just the DC component ∫B₀ dz — independent of R/σ (Boccaletti 1970;
Domcke, Garcia-Cely & Lee 2025, arXiv:2507.16609, Section 4.2).

The thin-lens / thick-lens distinction applies when Δk ≠ 0 (plasma detuning, axion mass),
creating a coherence length L_coh = 1/|Δk| that limits coherent accumulation. Without plasma,
L_coh → ∞ and the formula holds for all R (see [gertsenshtein_localized.md](gertsenshtein_localized.md)).

See [gertsenshtein_localized.md](gertsenshtein_localized.md) for the full physics derivation,
regime conditions, terminology mapping to literature, and open questions.

**Scripts**:
- `examples/gertsenshtein/run_localized.sh` — single run
- `examples/gertsenshtein/sweep_profile.sh` — 2D sweep over Bpeak × R

**Note**: Background B₀(z) is externally imposed — see [background_fields.md](background_fields.md)
for the validity discussion.

**References**: Boccaletti et al. (1970, Nuovo Cimento 70B:129); Raffelt & Stodolsky (1988,
PRD 37:1237); Domcke, Garcia-Cely & Lee (2025, arXiv:2507.16609)

### Target 5: Magnetar/FRB Application (Phase F3)

- Dipolar B(r) ∝ 1/r³ in 1+1D radial setup (spherical coordinates)
- Domain [R_NS, r_max] — avoids r=0 singularity
- Inner BC: Robin (impedance-matched) for graviton, Dirichlet for photon (conducting surface)
- Inward-propagating graviton wavepacket from large r
- Compare conversion efficiency with McDonald & Ellis (2024, arXiv:2406.18634)
- **Reference**: Kushwaha et al. (2022, arXiv:2202.00032); Domcke, Garcia-Cely & Lee (2025, arXiv:2507.16609)

## Implementation in TIDAL

The Gertsenshtein example is implemented using the **multi-field perturbation** pipeline extension (see [multi_field_perturbation.md](multi_field_perturbation.md)).

### TOML Configuration

See `examples/gertsenshtein/theory.toml` for the full configuration:
- Einstein-Maxwell Lagrangian with derived field F
- Metric perturbation h (TT gauge) + EM perturbation a (Lorenz gauge)
- Background 4-potential Ā generating uniform B₀
- Plane-wave reduction to 1+1D

### Pipeline Flow

```
theory.toml → _derive.py generates .wls with:
  1. DefTensorPerturbation for A → Ā + εa
  2. SetupMetricPerturbation for g → η + εh
  3. MakeRule expansion: F → dA (before perturbation)
  4. Perturbation[√|g| L, 2] + ExpandPerturbation
  5. Drop LI[2] terms, replace LI[1] → field tensors
  6. xCoba evaluates background: CD[-a][A[-b]] → F̄ = -B₀
  7. VarD for each dynamical field → coupled EOM
  8. DecomposeToComponents → JSON spec
→ JSON spec → simulate → measure conversion
```

## Future Extensions

1. **Plasma detuning** (Phase F1): ~~Proca-type mass term~~ — blocked by gauge-potential coupling artifact (see Target 3 above). Requires gauge-invariant mass mechanism for the transverse photon sector.
2. **Localized B-field scattering** (Phase F2): ✓ Implemented — `theory_localized.toml`, `run_localized.sh`, `sweep_profile.sh`. Validates Boccaletti formula P = sin²(κ × Bpeak × R × √(π/2)).
3. **Magnetar/FRB application** (Phase F3): Dipolar B(r), 1+1D radial setup, Gertsenshtein-Zel'dovich mechanism
4. **Non-minimal couplings** (Phase F5): Chern-Simons gravity (parity-violating P_L ≠ P_R), ξRF² curvature-EM coupling. Ref: Kushwaha & Jain (2024, arXiv:2410.07338)
5. **Axion-photon-graviton mixing** (future): Three-way mixing matrix with `g_{aγγ} a F F̃`. Ref: Dandoy, Lella et al. (2024)
6. **Photon-torsion conversion** (ultimate goal): Obukhov et al. (2024, arXiv:2410.01355) — axial torsion waves coupled to EM. Separate branch after Gertsenshtein program matures.
