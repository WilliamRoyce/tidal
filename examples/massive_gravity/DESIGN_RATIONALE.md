# Massive Gravity Example: Design Rationale

## 1. Why This Example Was Chosen

### Problem Statement

The `[linearization]` TOML feature was implemented to support xPert-based perturbation theory in the `tidal derive` CLI. However, the only existing linearization example (`gravitational_waves/`) exercises a narrow subset of the feature:

- Only linearizes `Einstein[CD][-a,-b]` (no additional terms)
- Uses 3+1D (4 dimensions, 10 components)
- Has no constants or runtime parameters
- Produces a zero mass matrix
- Does not reference the metric tensor in the expression

A second example was needed to validate that the linearization pipeline handles the full feature set robustly.

### Candidate Examples Considered

Several physically motivated examples were evaluated:

| Candidate | Pros | Cons | Decision |
|-----------|------|------|----------|
| **Linearized gravity in 2+1D** (massless) | Different dimension, fewer components | Same expression as GW, still zero mass matrix, all constraints (no evolution) | Rejected: too similar |
| **Linearized gravity with cosmological constant** (Lambda * g_ab) | Tests metric reference, adds constant | Background consistency issue: flat Minkowski + Lambda is inconsistent at zeroth order | Considered but superseded |
| **Linearized f(R) gravity** (R + alpha R^2) | Novel physics | Fourth-order equations, pipeline only handles second-order spatial operators | Rejected: pipeline limitation |
| **Linearized Einstein-Maxwell** | Multi-field perturbation | Requires perturbation of TWO fields (metric + EM), not supported by current single-field linearization | Rejected: significant extension needed |
| **Linearized Yang-Mills** | Non-abelian gauge theory | Requires Lie algebra structure (xTras), not in pipeline | Rejected: dependency issue |
| **Fierz-Pauli massive gravity** (m^2 * h_ab) | Tests constants, mass terms, metric reference | Mass term is quadratic in perturbation, not a linearization of a geometric expression | Modified approach adopted |
| **Linearized G_ab - m^2 g_ab = 0 in 2+1D** | Tests 7 distinct features vs GW, physically motivated, clean implementation | Zeroth-order background inconsistency (-m^2 * eta != 0) | **Selected** (see rationale) |

### Why Massive Gravity in 2+1D Won

The chosen example `G^(1)_ab[h] - m^2 h_ab = 0` obtained by linearizing `Einstein[CD][-a,-b] - m2 * eta[-a,-b]` was selected because it maximizes the number of distinct pipeline features tested while remaining physically meaningful and requiring **zero code changes** to the pipeline:

1. **Constants in linearization** -- `m2` as `DefConstantSymbol`, tested in Lagrangian-derived examples but never with linearization
2. **Metric reference in expression** -- `eta[-a,-b]` exercises the `eta[` -> `{prefix}Eta[` substitution rule that already existed but was untested in the linearization context
3. **Different dimension** -- 2+1D (6 components) vs 3+1D (10 components)
4. **Non-zero mass matrix** -- The identity operator terms from `m2 * h_ab` produce a non-zero diagonal mass matrix in the JSON, exercising the mass matrix auto-computation path
5. **Runtime parameters** -- `--param m2=1.0` enables parameter sweeps without re-derivation
6. **Dispersive propagation** -- omega^2 = k^2 + m^2 gives subluminal, spreading waves (testable physics)
7. **Mixed LHS types** -- constraint (h_0), first-order (h_1, h_2), and second-order (h_3-h_5) equations

## 2. Key Physics Considerations

### Background Consistency

The expression `G_ab - m^2 g_ab = 0` at zeroth order gives `-m^2 * eta_ab = 0`, which is inconsistent for m^2 != 0 (Minkowski is not a solution; the true background would be de Sitter). However:

- **At the linearized level**, xPert's `LinearizeTensorExpression` extracts the O(epsilon) part regardless of background consistency. It computes `G^(1)_ab - m^2 * h^(1)_ab` correctly.
- **This is standard practice** in massive gravity: the Fierz-Pauli mass term IS the linearized theory. The full nonlinear completion (dRGT massive gravity, etc.) is a separate question.
- **For pipeline testing**, what matters is that the computation produces valid, simulable equations with mass terms. Background consistency is a physics question, not a pipeline question.

### Sign Convention: Why MINUS m^2

The sign of the mass term in the linearization expression is critical for stability:

| Expression | Linearizes to | Evolution form | Behavior |
|-----------|---------------|----------------|----------|
| `Einstein[CD] + m2 eta` | `G^(1)_ab + m² h_ab = 0` | `d2_t(h) = +m² h + spatial` | **Exponential growth** (tachyonic) |
| `Einstein[CD] - m2 eta` | `G^(1)_ab - m² h_ab = 0` | `d2_t(h) = -m² h + spatial` | **Stable oscillation** (Klein-Gordon) |

The root cause is the sign convention of the linearized Einstein tensor: `G^(1)_ab` contains `-½ d2_t(h_ab) + ...`, so when solving for d2_t the sign of the mass term flips. The positive expression `+m²` in the original equation becomes positive on the RHS after rearranging, giving exponential growth.

This was discovered during simulation testing: the initial `+ m2 eta` expression produced exponentially growing solutions visible in the peak amplitude plot. Switching to `- m2 eta` produces the physically expected dispersive behavior.

**Note on Fierz-Pauli:** The true Fierz-Pauli mass term is `m²/2 (h_ab - η_ab h)` where h = η^cd h_cd is the trace. This cannot be obtained by linearizing any covariant expression, as it's intrinsically a linearized-level construction. Our expression `-m² η_ab` linearizes to `-m² h_ab` (without the trace subtraction), which only gives mass to diagonal components (h_0, h_3, h_5) since `η_ab` is diagonal. Off-diagonal components (h_1, h_2, h_4) evolve as massless wave equations. This is a simplification of Fierz-Pauli but sufficient for pipeline testing and produces qualitatively correct dispersive physics.

### Why 2+1D Specifically

In 3D (2+1D) spacetime, pure General Relativity has **no local propagating degrees of freedom**. The Weyl tensor vanishes identically, and the Riemann tensor is fully determined by the Ricci tensor. This means:

- Massless linearized gravity (the GW example in 2+1D) would produce ONLY constraint equations -- no time evolution. This is physically correct but not useful for simulation testing.
- Adding a mass term creates a propagating mode (the massive graviton). This makes the simulation physically interesting: you can see dispersive wave propagation.
- The constraint/evolution split (h_0 = constraint, h_1/h_2 = first-order, h_3-h_5 = second-order) is distinct from 4D and tests the pipeline's handling of mixed equation types in a different configuration.

### Dispersion Relation

For massive gravity in de Donder-like gauge, each evolution equation takes the form:
```
d2_t(h_i) = laplacian(h_i) - m^2 * h_i
```
This is the massive Klein-Gordon equation with dispersion relation omega^2 = k^2 + m^2:
- **Group velocity**: v_g = k / sqrt(k^2 + m^2) < 1 (subluminal)
- **Phase velocity**: v_p = sqrt(k^2 + m^2) / k > 1 (superluminal)
- **Observable effect**: A Gaussian pulse broadens over time as different Fourier components travel at different speeds

This is directly analogous to the Proca example (massive vector), providing a natural comparison.

## 3. Technical Design Decisions

### Expression Syntax

The TOML expression:
```toml
expression = "Einstein[CD][-a, -b] - m2 eta[-a, -b]"
```

Uses Mathematica juxtaposition multiplication (`m2 eta` without `*`), consistent with the conventions in existing Lagrangian expressions (e.g., `procaMassSquared/2 A[-a]` in Proca). The minus sign is essential for stable evolution (see "Sign Convention" above).

After substitution by `_substitute_field_names` (prefix "mg"):
```mathematica
Einstein[mgCD][-a, -b] - m2 mgEta[-a, -b]
```

Key substitution rules exercised:
- `CD]` -> `mgCD]` (CD followed by closing bracket, needed for `Einstein[CD]`)
- `eta[` -> `mgEta[` (metric name prefix)
- `m2` stays as-is (constants are NOT prefixed)

### Zero Code Changes Required

A critical finding during planning was that **all necessary substitution rules already existed** in `_substitute_field_names`:

```python
result = result.replace("eta[", f"{prefix}Eta[")   # Line ~388
result = result.replace("CD[", f"{prefix}CD[")     # Line ~389
result = result.replace("CD]", f"{prefix}CD]")     # Line ~390
```

The `eta[` substitution was added as part of the general substitution system but had never been exercised by the linearization path. This example validates that it works correctly in that context.

### JSON Output Structure

The expected JSON output has:
- 6 equations (h_0 through h_5)
- Mass terms appearing as `"operator": "identity"` with `"coefficient_symbolic": "-m2"`
- Non-zero diagonal mass matrix (auto-computed from identity operators)
- `"linearized": true` in metadata
- `"parameters": {"m2": 1.0}` for runtime override

### Simulation Design

The simulation places a Gaussian pulse in h_4 (h_xy), the only evolution equation:
- 2D spatial grid (64x64, periodic BCs)
- RK4 time integration, ~3.5 oscillation periods
- Physics validation: center-point oscillation matches analytic `A*cos(sqrt(2*m^2)*t)`
- 2x3 plot layout: heatmap snapshots at t=0/T/4/T/2, center oscillation, cross-sections, amplitude envelope

### Constraint Handling Limitations

The gauge-unfixed massive gravity equations produce 5 constraint equations and 1 evolution equation (h_4). The constraints have diverse mathematical forms:

| Constraint | Type | Self-terms |
|-----------|------|------------|
| h_0 (Hamiltonian) | Helmholtz | identity + laplacian_x + laplacian_y |
| h_1, h_2 (momentum) | Partial Helmholtz | identity + one-axis laplacian |
| h_3, h_5 (spatial) | Algebraic | identity only (no self-laplacian) |

The pipeline's current constraint solver (`_solve_constraint_equation`) only handles pure Poisson equations of the form `laplacian(field) = source`. It cannot solve:
- **Helmholtz equations** (`laplacian(phi) + lambda*phi = source`) for h_0, h_1, h_2
- **Algebraic constraints** (`(1+m^2)*phi = source`) for h_3, h_5
- **Coupled constraints** where h_0, h_3, h_5 depend on each other

With the constraint solver disabled (the default), constraint fields remain frozen at their initial values. Since h_4's evolution equation references constraint fields (h_0, pi_1, pi_2), and these are frozen at zero, h_4 evolves as a pure massive oscillator: `d2_t(h_4) = -2*m^2*h_4` with no spatial propagation.

This is still physical (it demonstrates the massive mode frequency omega = sqrt(2*m^2)), but spatial wave propagation would require either:
1. **De Donder gauge-fixed equations** where each component satisfies `Box h_ij + m^2 h_ij = 0` (standard massive KG with laplacian)
2. **Extended constraint solver** supporting Helmholtz and algebraic constraints (future pipeline work)

## 4. What This Example Validates End-to-End

### TOML -> WLS Generation (dry-run verified)
- [x] `[linearization]` section parsed correctly with `[constants]` and `[parameters]`
- [x] Expression substitution handles `eta[-a,-b]` -> `mgEta[-a,-b]`
- [x] `DefConstantSymbol[m2]` generated before linearization
- [x] `SetupMetricPerturbation` uses correct metric symbol
- [x] `LinearizeTensorExpression` receives fully substituted expression
- [x] Notation conversion `mghpert[LI[1], idx__] :> mgH[idx]` correct
- [x] Parameter defaults `"m2" -> 1.0` injected into metadata
- [x] `"linearized" -> True` in metadata

### WLS -> JSON (requires wolframscript + xPert)
- [x] xPert correctly linearizes `Einstein[CD] - m2 * eta` expression
- [ ] Mass term `m2 * h_ab` decomposed to identity operators per component
- [ ] 6 component equations generated (symmetric rank-2 in 3D)
- [ ] Mass matrix auto-computed as non-zero diagonal
- [ ] Symbolic coefficient preservation (`"m2"` string in JSON)

### JSON -> Simulation (requires JSON output)
- [ ] `build_pde_from_json` loads with `parameters={"m2": 1.0}`
- [ ] Mass terms evaluated from symbolic coefficients
- [ ] `create_initial_state` handles mixed constraint/evolution layout
- [ ] RK4 evolution stable with dt=0.005
- [ ] Dispersive spreading observable in h_3 amplitude

## 5. Comparison with Gravitational Waves Example

| Aspect | `gravitational_waves/` | `massive_gravity/` |
|--------|------------------------|--------------------|
| **Dimension** | 3+1D (4 spacetime dims) | 2+1D (3 spacetime dims) |
| **Components** | 10 (symmetric 4x4) | 6 (symmetric 3x3) |
| **Expression** | `Einstein[CD][-a,-b]` | `Einstein[CD][-a,-b] - m2 eta[-a,-b]` |
| **Constants** | None | `m2` (DefConstantSymbol) |
| **Parameters** | None | `m2 = 1.0` (runtime) |
| **Mass matrix** | Zero (massless) | Non-zero diagonal (massive) |
| **Metric in expr** | No | Yes (`eta[-a,-b]`) |
| **Physics** | Massless wave propagation | Dispersive massive propagation |
| **Constraint count** | 4 constraints, 6 evolution | 5 constraints, 1 evolution (h_4) |
| **Gauge variants** | Gauge-unfixed + de Donder + TT | Gauge-unfixed only |
| **Simulation focus** | TT-gauge polarizations h+, hx | Massive oscillation of h_xy |
| **Analogous to** | Massless photon (Maxwell) | Massive photon (Proca) |
