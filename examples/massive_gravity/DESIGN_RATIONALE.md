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
| **Linearized G_ab + m^2 g_ab = 0 in 2+1D** | Tests 7 distinct features vs GW, physically motivated, clean implementation | Zeroth-order background inconsistency (m^2 * eta != 0) | **Selected** (see rationale) |

### Why Massive Gravity in 2+1D Won

The chosen example `G^(1)_ab[h] + m^2 h_ab = 0` obtained by linearizing `Einstein[CD][-a,-b] + m2 * eta[-a,-b]` was selected because it maximizes the number of distinct pipeline features tested while remaining physically meaningful and requiring **zero code changes** to the pipeline:

1. **Constants in linearization** -- `m2` as `DefConstantSymbol`, tested in Lagrangian-derived examples but never with linearization
2. **Metric reference in expression** -- `eta[-a,-b]` exercises the `eta[` -> `{prefix}Eta[` substitution rule that already existed but was untested in the linearization context
3. **Different dimension** -- 2+1D (6 components) vs 3+1D (10 components)
4. **Non-zero mass matrix** -- The identity operator terms from `m2 * h_ab` produce a non-zero diagonal mass matrix in the JSON, exercising the mass matrix auto-computation path
5. **Runtime parameters** -- `--param m2=1.0` enables parameter sweeps without re-derivation
6. **Dispersive propagation** -- omega^2 = k^2 + m^2 gives subluminal, spreading waves (testable physics)
7. **Mixed LHS types** -- constraint (h_0), first-order (h_1, h_2), and second-order (h_3-h_5) equations

## 2. Key Physics Considerations

### Background Consistency

The expression `G_ab + m^2 g_ab = 0` at zeroth order gives `m^2 * eta_ab = 0`, which is inconsistent for m^2 != 0. This means Minkowski spacetime is not a solution of the full nonlinear equations. However:

- **At the linearized level**, xPert's `LinearizeTensorExpression` extracts the O(epsilon) part regardless of background consistency. It computes `G^(1)_ab + m^2 * h^(1)_ab` correctly.
- **This is standard practice** in massive gravity: the Fierz-Pauli mass term IS the linearized theory. The full nonlinear completion (dRGT massive gravity, etc.) is a separate question.
- **For pipeline testing**, what matters is that the computation produces valid, simulable equations with mass terms. Background consistency is a physics question, not a pipeline question.

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
expression = "Einstein[CD][-a, -b] + m2 eta[-a, -b]"
```

Uses Mathematica juxtaposition multiplication (`m2 eta` without `*`), consistent with the conventions in existing Lagrangian expressions (e.g., `procaMassSquared/2 A[-a]` in Proca).

After substitution by `_substitute_field_names` (prefix "mg"):
```mathematica
Einstein[mgCD][-a, -b] + m2 mgEta[-a, -b]
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
- Mass terms appearing as `"operator": "identity"` with `"coefficient_symbolic": "m2"` or `"-m2"`
- Non-zero diagonal mass matrix (auto-computed from identity operators)
- `"linearized": true` in metadata
- `"parameters": {"m2": 1.0}` for runtime override

### Simulation Design

The simulation follows the Proca simulation pattern:
- Gaussian pulse in h_3 (h_xx, a spatial-spatial component)
- 2D spatial grid (64x64, periodic BCs)
- RK4 time integration
- Physics validation: peak amplitude decreases over time (dispersive spreading)
- 2x2 plot layout: initial heatmap, final heatmap, amplitude decay, cross-section comparison

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
- [ ] xPert correctly linearizes `Einstein[CD] + m2 * eta` expression
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
| **Expression** | `Einstein[CD][-a,-b]` | `Einstein[CD][-a,-b] + m2 eta[-a,-b]` |
| **Constants** | None | `m2` (DefConstantSymbol) |
| **Parameters** | None | `m2 = 1.0` (runtime) |
| **Mass matrix** | Zero (massless) | Non-zero diagonal (massive) |
| **Metric in expr** | No | Yes (`eta[-a,-b]`) |
| **Physics** | Massless wave propagation | Dispersive massive propagation |
| **Constraint count** | 4 constraints, 6 evolution | 1 constraint, 2 first-order, 3 evolution |
| **Gauge variants** | Gauge-unfixed + de Donder + TT | Gauge-unfixed only |
| **Simulation focus** | TT-gauge polarizations h+, hx | Dispersive spreading of h_xx |
| **Analogous to** | Massless photon (Maxwell) | Massive photon (Proca) |
