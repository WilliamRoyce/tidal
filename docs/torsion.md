# Torsion in Poincaré Gauge Theory

## Overview

TIDAL supports Poincaré gauge theory (PGT) Lagrangians with propagating torsion. The torsion tensor T^a_{bc} is the field strength of the translational gauge symmetry, introduced via xAct's `DefCovD` with `Torsion -> True`. Torsion components are perturbed alongside the metric perturbation using `DefTensorPerturbation`, enabling the full graviton-torsion mixing problem.

## PGT Framework

In Poincaré gauge theory, gravity is a gauge theory of the Poincaré group with two independent gauge fields:

- **Tetrad** (vierbein) e^a_μ: gauges translations → encodes the metric via g_μν = η_{ab} e^a_μ e^b_ν
- **Spin connection** ω^{ab}_μ: gauges Lorentz rotations

The corresponding field strengths are:

- **Torsion**: T^a_{μν} = ∂_μ e^a_ν - ∂_ν e^a_μ + ω^a_{bμ} e^b_ν - ω^a_{bν} e^b_μ
- **Curvature**: R^{ab}_{μν} (standard Riemann-Cartan curvature)

The Riemann-Cartan Ricci scalar R̃ includes torsion contributions and is related to the Levi-Civita Ricci scalar R by the identity:

```
R̃ = R - 1/4 T_{abc}T^{abc} - 1/2 T_{abc}T^{bac} + T^a_{ab}T^{cb}_c + total derivative
```

(Shapiro 2002, hep-th/0103093, eq. 2.17; Hehl et al. 1976)

## TOML Configuration

### Enabling Torsion

```toml
[spacetime]
dimension = 4
metric = "minkowski"
torsion = true   # Creates CDT with DefCovD[..., Torsion -> True]
```

This auto-creates:
- `TorsionCDT[a, -b, -c]`: torsion tensor (antisymmetric in -b, -c)
- `RicciScalarCDT[]`: full Riemann-Cartan Ricci scalar
- Other curvature tensors with torsion contributions

### Writing the Lagrangian

The Lagrangian must be written in **decomposed form** using `RicciScalarCD[]` (Levi-Civita, NOT `RicciScalarCDT[]`) plus explicit `TorsionCDT[...]` terms:

```toml
[lagrangian]
expression = "(1/kappa^2) RicciScalarCD[] + alpha1 TorsionCDT[-a,-b,-c] eta[a,d] eta[b,e] eta[c,f] TorsionCDT[-d,-e,-f] + alpha2 TorsionCDT[-a,-b,-c] eta[b,d] eta[a,e] eta[c,f] TorsionCDT[-d,-e,-f] + alpha3 TorsionCDT[b,-b,-a] eta[a,d] TorsionCDT[c,-c,-d]"
```

**Why decomposed form?** xPert's built-in perturbation formulas assume the Levi-Civita connection. Using `RicciScalarCDT[]` directly would require decomposing it via `ChangeCurvature`, but `ChangeTorsion` converts torsion → contortion (a derived quantity that `DefTensorPerturbation` cannot perturb). The decomposed form avoids this issue entirely.

The three T² invariants correspond to the three irreducible pieces of torsion under SO(1,3):
- `T_{abc}T^{abc}`: couples to all three irreducible sectors
- `T_{abc}T^{bac}`: the "mixed" contraction
- `T^a_{ab}T^{cb}_c`: the trace sector

### Linearization

```toml
[linearization]
perturbation_field = "h"   # Metric perturbation (standard xPert)
```

Torsion perturbation is **auto-registered** when `torsion = true`. No `[[linearization.matter_perturbations]]` entry is needed — the pipeline automatically calls `DefTensorPerturbation` on `TorsionCDT` with background = 0.

## Torsion Irreducible Decomposition

Torsion T_{abc} (24 components in 4D) decomposes into three irreducible pieces under SO(1,3):

| Piece | Name | Components | Description |
|-------|------|-----------|-------------|
| ^(1)T_{abc} | Tensor (tentor) | 16 | Traceless, tensor part |
| ^(2)T_a = T^b_{ba} | Trace (trator) | 4 | Vector trace |
| ^(3)T_a = ε_{abcd}T^{bcd} | Axial (axitor) | 4 | Totally antisymmetric → pseudoscalar |

The PGT Lagrangian parameters α₁, α₂, α₃ control the mass and kinetic terms for each sector independently.

### Propagating Modes

| PGT Lagrangian | Propagating torsion DOF |
|---------------|------------------------|
| R only (Einstein-Cartan) | 0 — torsion is non-propagating |
| R + α₁T² | Massive tensor + other modes |
| R + α₃T²_trace | 1 pseudoscalar (axial) |
| R + α₁T² + α₂T² + α₃T² | Up to 5 modes |

**Key result**: Pure Einstein-Cartan (L = R with torsion, no T² terms) gives NO torsion-Gertsenshtein effect — torsion is non-propagating and vanishes for spinless EM (theorem, not approximation). The T² kinetic terms are essential for propagating torsion.

## Implementation Architecture

### What Happens During Derivation

1. `torsion = true` → `DefCovD[CDT, Torsion -> True, FromMetric -> metric]`
2. Lagrangian parsed with `RicciScalarCD[]` + `TorsionCDT[...]` terms
3. `SetupMetricPerturbation` for metric h
4. Auto `DefTensorPerturbation` for `TorsionCDT` (background = 0)
5. xPert: L^(2) = second-order perturbation in both h and TorsionCDT
6. VarD: derive EOM for h components and TorsionCDT components
7. Component decomposition → plane-wave reduction → constraint elimination

### Why ChangeCurvature/ChangeTorsion Were NOT Used

Investigation revealed:
- `ChangeCurvature[L, CDT, CD]` correctly decomposes R̃ → R + contortion terms
- But the contortion terms involve `Christoffel[CD,CDT]` — a derived xAct object
- `ChangeTorsion[L, CDT, CD]` converts TorsionCDT → Christoffel[CD,CDT] (the WRONG direction)
- `DefTensorPerturbation` requires fundamental tensors, not derived Christoffel differences
- xAct does not provide the inverse conversion (contortion → torsion) natively

The decomposed-form Lagrangian avoids this issue entirely.

## References

- Blagojevic, M. & Hehl, F.W. (2013). *Gauge Theories of Gravitation*. Imperial College Press.
- Hehl, F.W. et al. (1976). "General Relativity with Spin and Torsion." *Rev. Mod. Phys.* 48:393. DOI: [10.1103/RevModPhys.48.393](https://doi.org/10.1103/RevModPhys.48.393)
- Shapiro, I.L. (2002). "Physical Aspects of the Space-Time Torsion." *Phys. Rept.* 357:113. arXiv: [hep-th/0103093](https://arxiv.org/abs/hep-th/0103093)
- Hayashi, K. & Shirafuji, T. (1979). "New General Relativity." *Phys. Rev. D* 19:3524.
- Nikiforova, V. et al. (2009). "Stability of the Massive Torsion Modes." arXiv: [0905.4007](https://arxiv.org/abs/0905.4007)
- Barker, W.E.V. (2023). "HiGGS: Hamiltonian in Gauge Gravity Solver." *Eur. Phys. J. C* 83:228. arXiv: [2206.00658](https://arxiv.org/abs/2206.00658)
- xAct: [DefCovD with Torsion](http://www.xact.es/Documentation/HTML/HTMLLinks/xTensor/DefCovD.nb.html)
