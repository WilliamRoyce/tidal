# Coupled Proca: Two Massive Vector Fields in a Cavity

## Purpose

This example is a **constraint solver stress test**, designed to exercise
non-periodic code paths that are not covered by other examples:

| Feature | Coverage before | This example |
|---------|----------------|--------------|
| Coupled FFT solve (periodic) | massive_gravity | - |
| Matrix/sparse solver (non-periodic) | electrostatics (uncoupled) | **Coupled Helmholtz** |
| Gauss-Seidel coupled iteration | 1 unit test (algebraic only) | **Real spatial operators** |
| Dirichlet BCs on constraints | electrostatics (uncoupled) | **Coupled system** |
| Two vector fields in pipeline | None | **First example** |

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

Non-periodic Dirichlet (rectangular cavity, fields vanish on walls):
- Forces the **sparse matrix solver** path (not FFT)
- Forces **Gauss-Seidel iteration** for coupled constraints
- Physical interpretation: two vector fields confined in a box

## Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| mA^2 | 1.0 | Helmholtz scale for A |
| mB^2 | 2.0 | Different Helmholtz scale for B |
| g | 0.5 | Cross-coupling strength |

Different masses test heterogeneous constraint operators.

## Running

```bash
# Simulation script
uv run python examples/coupled_proca/simulation.py

# Or via CLI
./examples/coupled_proca/run.sh

# Re-derive from Lagrangian (requires wolframscript)
tidal derive examples/coupled_proca/theory.toml
```

## Validation

1. **Constraints activate**: A_0, B_0 become non-zero when momenta are non-zero
2. **Coupling works**: B_1, B_2 become non-zero (energy transfer A -> B)
3. **Stability**: All fields remain finite and bounded
4. **Gauss-Seidel converges**: No convergence warnings
5. **Dirichlet respected**: Boundary values stay near zero
