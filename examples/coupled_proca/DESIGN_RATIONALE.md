# Coupled Proca: Two Massive Vector Fields in 2+1D

## Purpose

This example is a **coupled constraint solver test**, designed to exercise
multi-field constraint solving with two distinct Helmholtz scales:

| Feature | Coverage before | This example |
|---------|----------------|--------------|
| Coupled FFT solve (periodic) | massive_gravity (single scale) | **Two different Helmholtz scales** |
| Cross-field identity coupling | scalar_vector_coupling | **Two vector fields** |
| Two vector fields in pipeline | None | **First example** |
| Measurement module on vectors | scalar_vector_coupling | **Group conversion A→B** |

## Lagrangian

Two massive vector fields A and B with cross-coupling:

```
L = -1/4 F^A_{ab} F^{A,ab} - 1/4 F^B_{ab} F^{B,ab}
    - mA^2/2 A_a A^a - mB^2/2 B_a B^a + g A_a B^a
```

## Component Structure (6 fields in 2+1D)

| Field | Time order | Type | Operators |
|-------|-----------|------|-----------|
| A_0 | 0 | Constraint (Helmholtz) | laplacian - mA^2 identity + g identity(B_0) |
| A_1 | 2 | Evolution | wave + mass + coupling |
| A_2 | 2 | Evolution | wave + mass + coupling |
| B_0 | 0 | Constraint (Helmholtz) | laplacian - mB^2 identity + g identity(A_0) |
| B_1 | 2 | Evolution | wave + mass + coupling |
| B_2 | 2 | Evolution | wave + mass + coupling |

The two constraints are **coupled** through the g identity cross-terms.

## Boundary Conditions

Periodic on both axes:
- Enables the **coupled FFT solver** path for constraints
- Energy measurement drift converges as O(dx^2) (~3% at 20x20, <1% at 48x48)
- Physical interpretation: infinite-domain physics with periodic images

## Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| mA^2 | 1.0 | Helmholtz scale for A |
| mB^2 | 2.0 | Different Helmholtz scale for B |
| g | 0.5 | Cross-coupling strength |

Different masses test heterogeneous constraint operators.

## Running

```bash
# Full pipeline (derive + simulate + measure + plot)
cd examples/coupled_proca && bash run.sh

# Or re-derive from Lagrangian only (requires wolframscript)
tidal derive examples/coupled_proca/theory.toml
```

## Measurement

`tidal measure` computes:
- **Group conversion**: `A_1 -> {A_2, B_1, B_2}` energy transfer via `compute_group_conversion`
- **Spectral conversion** P(k,t) and **dispersion** omega(k)
- **Energy conservation diagnostics**: O(dx^2) convergent with periodic BCs

## Validation

1. **Constraints activate**: A_0, B_0 become non-zero when momenta are non-zero
2. **Coupling works**: B_1, B_2 become non-zero (energy transfer A -> B)
3. **Stability**: All fields remain finite and bounded
4. **FFT solver converges**: Coupled constraints solved without warnings
5. **Energy drift convergent**: O(dx^2) with resolution, ~3% at 20x20
