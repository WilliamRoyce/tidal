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

**Note on P&R discrepancy:** Palessandro & Rothman (2023, Eq. 26) quote
P = sin²(√G · B₀ · D) = sin²(κB₀D/(4√π)), which differs from our result by a
factor of 2√π ≈ 3.54. The source of this discrepancy is under investigation. Our
formula is derived directly from the E-L equations of (1/κ²)R - (1/4)F² with
g = η + h and verified numerically.

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
| [**Palessandro & Rothman, PDU 40:101187**](https://arxiv.org/abs/2301.02072) | **2023** | **Simple derivation from Lagrangian — our primary validation target** |
| [Palessandro, arXiv:2405.01407](https://arxiv.org/abs/2405.01407) | 2024 | Most complete reference: full Lagrangian, QED corrections, mixing matrix |

### Extensions and Modern Applications

| Ref | Year | Key Contribution |
| --- | ---- | ---------------- |
| [Domcke & Garcia-Cely, JCAP 05:051](https://arxiv.org/abs/2312.17636) | 2024 | Inverse Gertsenshtein as HFGW probe |
| [Dandoy et al., arXiv:2406.17853](https://arxiv.org/abs/2406.17853) | 2024 | Full 4-component mixing with plasma/QED effects |
| [Obukhov et al., arXiv:2410.01355](https://arxiv.org/abs/2410.01355) | 2024 | Photon-torsion wave conversion in Poincaré gauge theory |
| [arXiv:2507.02362](https://arxiv.org/abs/2507.02362) | 2025 | EM-torsion coupling near black holes |

### Numerical Methods

| Ref | Year | Key Contribution |
| --- | ---- | ---------------- |
| [Berlin et al., arXiv:2405.08865](https://arxiv.org/abs/2405.08865) | 2024 | Numerical analysis of resonant axion-photon mixing — directly analogous methodology |

## Validation Targets

All equations are **derived from the Lagrangian** by the TIDAL pipeline — listed here only for comparison.

### Target 1: Vacuum Conversion Probability

**Source**: Palessandro & Rothman (2023), Eq. 26.

- Uniform B₀, vacuum (no plasma), massless graviton
- Plane wave at known wavenumber k propagating along z
- Measure P(t) via `tidal measure --what conversion`
- Compare: `P = sin²(κ B₀ D / 2)`
- **Acceptance**: Agreement to within 1%

### Target 2: B₀ Sweep

- Sweep B₀ from weak to strong coupling
- Plot P_max vs B₀
- Compare against analytical curve
- Check oscillation length scaling: L_osc ∝ 1/B₀

### Target 3: Detuned Conversion (Future — Phase F1)

- Add effective plasma mass m_γ to EM perturbation
- Sweep ω_p at fixed B₀
- Compare with full detuned formula
- Test MSW-like resonance when ω_p² = m_g²
- **Acceptance**: Agreement to within 5%

### Target 4: Finite-Region Scattering (Phase E)

- Localized B₀(z) with Gaussian profile
- Wave packet enters, converts, exits
- Sweep region width R
- Compare against Boccaletti et al. integral formula

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

1. **Plasma/QED effects** (Phase F1): Effective photon mass from plasma frequency
2. **Chern-Simons gravity** (Phase F2): Parity-violating h-a coupling, different L/R conversion rates
3. **Axion-photon-graviton mixing** (Phase F3): Three-way mixing matrix with `g_{aγγ} a F F̃`
4. **Photon-torsion conversion** (Phase F4): Obukhov et al. (2024) — axial torsion waves coupled to EM via Chern-Simons-like term. This is the ultimate project goal.
5. **Non-abelian Yang-Mills** (Phase F5): Autocatalytic graviton production without external B field (Palessandro & Rothman 2023)
