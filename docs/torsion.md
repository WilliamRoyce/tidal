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

Add a `[torsion]` section to the TOML config. This is an optional, standalone section — existing theories without torsion are unaffected.

```toml
[spacetime]
dimension = 4
metric = "minkowski"

[torsion]
perturbation_name = "t"     # name for the linearized torsion field in equations/JSON
# background = "Tbar"       # future: non-zero background torsion
# irreducible = "axial"     # future: restrict to specific torsion sector
```

The `[torsion]` section:
1. Extends the connection: defines CDT with `DefCovD[..., Torsion -> True]`
2. Auto-creates `TorsionCDT[a, -b, -c]` (antisymmetric in -b, -c) for use in the Lagrangian
3. Registers a perturbation field named by `perturbation_name` (e.g., `t_0`, `t_1`, ..., `t_8` in the JSON)
4. Handles background torsion zeroing (flat Minkowski) and xPert label management

**Architecture**: The `[torsion]` section is separate from `[spacetime]` because in PGT, the metric (from tetrad) and torsion (from connection) are independent gauge fields. The `[spacetime]` section defines the metric and coordinates; `[torsion]` extends the connection.

**Reserved names**: The `perturbation_name` must not collide with field names, constant names, or built-in operators (CD, CDT, eta, Ricci, Torsion, etc.). The validation system checks this automatically.

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

Torsion perturbation is **auto-registered** when a `[torsion]` section is present. No `[[linearization.matter_perturbations]]` entry is needed — the pipeline automatically:
1. Defines the perturbation field (e.g., `t[a,-b,-c]`) from `perturbation_name`
2. Calls `DefTensorPerturbation` to connect it to `TorsionCDT`
3. Handles the xPert label/truncation (see below)
4. Sets background torsion and contortion to zero (flat Minkowski)
5. Calls VarD w.r.t. the torsion perturbation field for the torsion EOM

### Perturbation Truncation Logic

xPert expands each field as: `Φ = Φ̄ + ε·δΦ^(1) + ε²·δΦ^(2) + ...`
- `LI[1]` = `δΦ^(1)` — the linear perturbation (dynamical field)
- `LI[2]` = `δΦ^(2)` — 2nd-order correction to the field itself (NOT a separate DOF)

The second-order Lagrangian L^(2) = ½δ²S contains three types of terms:
1. `(δΦ^(1))²` — quadratic products of 1st-order perturbations → **KEPT** (give the linear EOM)
2. `Φ̄ · δΦ^(2)` — background × 2nd-order correction → **ZERO** (background Φ̄ = 0)
3. `δΦ^(2)` alone — appears in L^(1) only, not L^(2)

Therefore `LI[2] → 0` is correct: the physical content is entirely in the `LI[1]` terms.
This applies identically to metric (h), torsion (t), and matter (a) perturbations.

### Background Contortion Zeroing

The contortion `K^a_{bc} = ½(T^a_{bc} + T_b^a_c - T_{bc}^a)` is algebraically
determined by torsion. For zero background torsion (T̄ = 0), the background
contortion K̄ = 0. In xPert's expansion:
- `Perturbation[K, 0]` = K̄ = 0 → zeroed via `ChristoffelCDCDT[__] :> 0`
- `Perturbation[K, 1]` = δK(t) → **preserved** (encodes the torsion wave)
- `Perturbation[K, 2]` → already dropped by `LI[2] → 0`

This is exactly analogous to zeroing `RicciScalarCD[]` (flat background curvature)
while keeping `Perturbation[RicciScalarCD[]]` (the graviton wave).

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

1. `[torsion]` section present → `DefCovD[CDT, Torsion -> True, FromMetric -> metric]`
2. Perturbation field defined: `DefTensor[{prefix}T[a,-b,-c], M, Antisymmetric[{-b,-c}]]`
3. `DefTensorPerturbation` connects `{prefix}tPert[LI[order],...]` to `TorsionCDT[...]`
4. Lagrangian parsed with `RicciScalarCD[]` + `TorsionCDT[...]` terms
5. `SetupMetricPerturbation` for metric h
6. xPert: L^(2) = second-order perturbation in both h and torsion
7. Truncation: `tPert[LI[2],...] → 0`, `tPert[LI[1],...] → T[...]`, `TorsionCDT[...] → 0`
8. VarD: derive EOM for both h and T components
9. Component decomposition → plane-wave reduction → constraint elimination

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
