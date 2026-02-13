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

## Measurement

`measure_coupling.py` uses the measurement module to compute:
- **Group conversion**: `A_1 -> {A_2, B_1, B_2}` energy transfer via `compute_group_conversion`
- **Per-component breakdown**: individual `A_1 -> A_2`, `A_1 -> B_1`, `A_1 -> B_2` probabilities
- **Energy conservation diagnostics**: with a relaxed threshold (0.5)

### Energy Conservation Caveat (Dirichlet + curl-curl)

The discrete `cross_derivative_xy` operator with Dirichlet ghost cells
is **not self-adjoint** at boundary cells:

```
M[(0,j), (0,j+1)] = +1/(4 dx dy)
M[(0,j+1), (0,j)] = -1/(4 dx dy)    ← opposite sign!
```

This makes the discrete curl-curl system non-Hamiltonian — no quadratic
energy functional is exactly conserved.  The ~30% energy drift occurs
even when using py-pde's own operators to compute the virial potential.
With **periodic BCs**, the same system conserves energy to machine
precision (~1e-10).

This is a fundamental discretization limitation of `cross_derivative`
operators on non-periodic grids, not a bug.  The measurement script
uses a relaxed threshold and documents this as an expected artifact.
Conversion probability measurements remain valid regardless.

## Validation

1. **Constraints activate**: A_0, B_0 become non-zero when momenta are non-zero
2. **Coupling works**: B_1, B_2 become non-zero (energy transfer A -> B)
3. **Stability**: All fields remain finite and bounded
4. **Gauss-Seidel converges**: No convergence warnings
5. **Dirichlet respected**: Boundary values stay near zero
