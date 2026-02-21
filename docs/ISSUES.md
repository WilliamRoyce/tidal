# Known Open Issues

## Non-Diagonal Kinetic Matrix in Linearized Gravity

**Status:** Open (blocking simulation of linearized gravity examples)
**Affects:** `massive_gravity`, `gravitational_waves` examples
**Priority:** High (unblocks Gertsenshtein effect simulation)
**Resolves with:** Phase B (Gauge Fixing) or explicit kinetic matrix inversion

### Problem

The canonical momentum pipeline computes field rates via the simple subtraction
`dq_i/dt = vel_i - pi_i`, which assumes a **unit (diagonal) kinetic coefficient**:
`pi_i = d_t q_i + spatial_terms`.

For linearized gravity, the quadratic Lagrangian `L^(2)` has cross-coupled time
derivatives of the form `(d_t h_{ab})(d_t h_{cd})` with various metric contractions.
This produces a **non-diagonal kinetic matrix** `K_{ij}`:

```
pi_i = K_{ij} * d_t q_j + spatial_terms(q)
```

When `K` is not the identity, the field rate `dq_i/dt = K^{-1}_{ij}(pi_j - spatial_j)`
requires matrix inversion. The current pipeline instead produces field_rates containing
`first_derivative_t` references to other fields, which the PDE builder rejects:

```
RuntimeError: Operator 'first_derivative_t' cannot be applied as a spatial operator.
```

### Root Cause

In the gauge-unfixed formulation of linearized Einstein equations, all 10 (or 6 in 2+1D)
components of `h_{ab}` appear as 2nd-order-in-time evolution equations. The kinetic
energy terms in `L^(2)_EH` couple different components' time derivatives:

```
L^(2) ~ (d_t h_xx)^2 + (d_t h_yy)^2 - (d_t h_xx)(d_t h_yy) + ...
```

This gives a kinetic matrix like:
```
K = [[1, 0, -1/2, ...],
     [0, 1,  0,   ...],
     [-1/2, 0, 1, ...],
     ...]
```

### Why Gauge Fixing Resolves This

In **de Donder (harmonic) gauge**, the linearized Einstein equations simplify to:

```
-1/2 Box h_ab + mass_terms = 0  (after trace-reversal)
```

where `Box = -d^2_t + Laplacian` is the flat-space d'Alembertian. This form has:
1. **Diagonal kinetic matrix** (each component has its own `d^2_t`)
2. **Explicit constraint structure** (h_0, h_0i become constraints)
3. **Standard wave equation form** that the PDE builder handles natively

### Resolution Options

1. **Phase B: Automatic gauge fixing** (recommended)
   - Implement de Donder/Lorenz gauge in the pipeline
   - Gauge fixing term `L_gf = -1/(2xi) (d^a h_ab - 1/2 d_b h)^2` diagonalizes kinetic matrix
   - Most physically correct approach; also exposes constraint structure
   - Already planned as Phase B in NEXT_PHASES.md

2. **Symbolic kinetic matrix inversion** (alternative)
   - Compute `K_{ij}` in Wolfram, invert symbolically, emit `K^{-1}` to JSON
   - More general but complex; may produce large expressions

3. **Runtime numeric matrix inversion** (fallback)
   - Build `K` numerically from JSON coefficients, invert with NumPy
   - Simple but may not handle position-dependent metrics

### Current State

The equations in `massive_gravity_3d.json` and `linearized_gravity.json` are
**physically correct** (contain full linearized Einstein tensor + mass terms from
the single-path L^(2) = Perturbation[L,2]/2 approach). The canonical Hamiltonian
is also correct (60 terms for massive gravity, 170 for gravitational waves).

Only the **field_rates** have the non-diagonal kinetic matrix issue, preventing
the PDE builder from constructing the time evolution. The equations, mass/coupling
matrices, and Hamiltonian are all valid.

### Related Issues

- Gauge-unfixed equations have no explicit constraint structure (all components
  are time_order=2). Gauge fixing will reveal Hamiltonian and momentum constraints.
- The `test_constraint_h0_nonzero` test in `test_constraint_solver.py` assumed
  h_0 was a constraint (from the old legacy linearization path). With the new
  Lagrangian-first approach, h_0 is correctly an evolution equation in the
  gauge-unfixed formulation.

### References

- Carroll, "Spacetime and Geometry" (2004), Ch. 7 (linearized gravity)
- Brizuela, Martin-Garcia, Mena Marugan (2009) (xPert perturbation theory)
- Fierz & Pauli (1939) (massive gravity mass term)
- Phase B plan in `docs/NEXT_PHASES.md`
