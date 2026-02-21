# Known Open Issues

## Non-Diagonal Kinetic Matrix in Linearized Gravity

**Status:** Resolved — K^{-1} symbolic inversion + Expand fix
**Affects:** Previously affected `massive_gravity`, `gravitational_waves` examples
**Priority:** Closed
**Resolved by:** Symbolic kinetic matrix inversion (K^{-1}) in canonical pipeline,
plus always-Expand fix in ExportJSON.wl

### Problem

The canonical momentum pipeline originally computed field rates via simple subtraction
`dq_i/dt = vel_i - pi_i`, which assumed a **unit (diagonal) kinetic coefficient**:
`pi_i = d_t q_i + spatial_terms`.

For linearized gravity, the quadratic Lagrangian `L^(2)` has cross-coupled time
derivatives of the form `(d_t h_{ab})(d_t h_{cd})` with various metric contractions.
This produces a **non-diagonal kinetic matrix** `K_{ij}`:

```
pi_i = K_{ij} * d_t q_j + spatial_terms(q)
```

When `K` is not the identity, the field rate `dq_i/dt = K^{-1}_{ij}(pi_j - spatial_j)`
requires matrix inversion. Without this, the pipeline produced field_rates containing
`first_derivative_t` references to other fields, which the PDE builder rejected.

### Resolution (Implemented)

**1. Symbolic K^{-1} inversion (gauge-independent)**

The canonical pipeline (`_wls_canonical_pipeline()` in `_derive.py`) now:
1. Computes the kinetic matrix `K_{ij} = D[pi_i, vel_j]` symbolically
2. Validates `det(K) != 0` (throws if singular)
3. Inverts symbolically via `Simplify[Inverse[K]]`
4. Emits field rates: `dq_i/dt = Sum_j K^{-1}_{ij} (pi_j - S_j)`
5. Verifies no residual `first_derivative_t` operators remain

This works for ALL gauge-unfixed theories without requiring any gauge choice:
- Proca (K = diag(-1, 1))
- Massive gravity (4x4 non-diagonal K)
- Gravitational waves (7x7 non-diagonal K, 3 constraint components)
- Elasticity (K = diag(rho, rho))

**2. Always-Expand fix in ExportJSON.wl**

Mathematica's `Total[]` in `EquationToJSONMultiField` could trigger auto-factoring,
collapsing separate linear terms into a single `Times` with multiple field heads.
`Expand` was only applied when `|lhsCoeff| != 1`, missing standard wave equations.

Fix: always call `Expand[rhs]` before `ParseMultiFieldRHS`, regardless of lhsCoeff.
Same defensive `Expand` added inside `ParseHamiltonianExpression`.

**3. EOM-based fast path for high-rank tensors**

For theories with raw component count > 30 (e.g., rank-3 antisymmetric in 4D:
4^3 = 64 raw components), `DecomposeScalarExpression` on the abstract Lagrangian
is prohibitively slow. The pipeline uses an EOM-based fast path that constructs
canonical structure directly from already-decomposed EOM (K=I assumed).

### Additional Options (for physics optimisation, not required for correctness)

1. **TT gauge** — reduces gravitational waves from 10 components to 2 physical
   polarisation modes (h_+, h_x). Available via `[[gauge]] type = "tt"`.

2. **De Donder (harmonic) gauge** — diagonalizes kinetic matrix and exposes
   constraint structure. Available via `[[gauge]] type = "de_donder"` (Phase B).

### Root Cause Analysis

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

### Current State

All 23 example JSON specs load, build PDE systems, and evolve correctly.
1016 Python tests pass, 0 failures, 0 xfails. The `_NON_EVOLVABLE_SPECS` set
in `test_example_jsons.py` is now empty.

### References

- Carroll, "Spacetime and Geometry" (2004), Ch. 7 (linearized gravity)
- Brizuela, Martin-Garcia, Mena Marugan (2009) (xPert perturbation theory)
- Fierz & Pauli (1939) (massive gravity mass term)
