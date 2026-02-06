# Chern-Simons Implementation Notes

## 📝 About This File

This is an **example-specific** implementation guide for the Chern-Simons 2+1D case.

**When to create similar files:**
- For any complex example requiring hybrid approaches (symbolic + manual)
- When implementing topological field theories with special tensor structures
- For examples with gauge symmetries or constraint equations
- When exploring new physics domains (Yang-Mills, BF theory, supergravity, etc.)

**Template for new example notes:**
1. Physics Background (Lagrangian, EOM, physical interpretation)
2. Implementation Status (what works, what doesn't, current limitations)
3. Wolfram Implementation Pattern (specific approaches for this example)
4. JSON Structure (how the equations are represented)
5. Python Simulation (grid setup, initial conditions, expected physics)
6. Future Automation Plan (steps to fully automate if currently hybrid)
7. References (scripts, JSON files, related examples)

**Integration with main docs:**
- Reference this file from `MEMORY.md` under "Example Implementations"
- Add error patterns to `troubleshooting.md` if they're general enough
- Link back to main docs for general patterns (this file is for specifics)

---

## Physics Background

Chern-Simons theory in 2+1D:
```
L = -1/4 F_μν F^μν + (κ/2) ε^μνρ A_μ ∂_ν A_ρ
```

Equations of motion:
```
∂_μ F^μν + κ ε^ναβ ∂_α A_β = 0
```

In Lorenz gauge (∂^μ A_μ = 0):
```
d²A_0/dt² = ∇²A_0 + κ(∂_x A_2 - ∂_y A_1)
d²A_1/dt² = ∇²A_1 + κ(∂_y A_0)  [time term -κ∂_t A_2 moves to LHS]
d²A_2/dt² = ∇²A_2 - κ(∂_x A_0)  [time term +κ∂_t A_1 moves to LHS]
```

## Implementation Status

### What Works (Fully Automated)
- ✅ 3D manifold setup (M3, 2+1D spacetime)
- ✅ Full Maxwell-Chern-Simons Lagrangian symbolic derivation
- ✅ **Automatic epsilon tensor evaluation** (EvaluateEpsilonComponents)
- ✅ Component decomposition for vector fields in 3D
- ✅ **Automatic gradient direction detection** (gradient_x vs gradient_y)
- ✅ **Cross-field reference detection** (A_0 → A_1 → A_2)
- ✅ Python simulation with gradient operators
- ✅ Energy transfer between components (CS coupling effect)

### Known Limitations
- ⚠️ Mixed time-space derivatives (d_t d_x) classified as "laplacian" in JSON
- ⚠️ Extra cross-Laplacian terms appear (from unsimplified gauge terms)
- ⚠️ Time derivative separation pattern needs 2+1D update for cleaner JSON

## Wolfram Implementation Pattern

### Fully Automated Approach (Current)
```mathematica
(* 1. Build full Lagrangian with epsilon tensor *)
MaxwellLagrangian = -1/2 CD3[-a][csA[-b]] eta3[a,c] eta3[b,d] CD3[-c][csA[-d]]
CSLagrangian = (kappa/2) * epsiloneta3[a, b, c] * csA[-a] * CD3[-b][csA[-c]]
FullLagrangian = MaxwellLagrangian + CSLagrangian

(* 2. Derive EOM - includes both Maxwell and CS contributions *)
eomFull = VarD[csA[-a], CD3][FullLagrangian]

(* 3. Decompose - EvaluateEpsilonComponents handles epsilon automatically *)
componentEqs = DecomposeToComponents[eomFull, csA[-a], cart3]

(* 4. Export - IdentifyGradientDirection detects x/y correctly *)
jsonStructure = BuildMultiFieldJSONStructure[fullEquations, metadata]
```

### Previous Hybrid Approach (No Longer Needed)
The old approach required manual CS term addition after decomposition.
This is no longer necessary with the automated epsilon evaluation.

### Why Automation Now Works
The epsilon tensor ε^μνρ creates complex index structures:
- xAct's `epsiloneta3[-a,-b,-c]` exists but component evaluation is non-trivial
- `ToBasis[chart][epsilon[...]]` produces symbolic basis components
- Need explicit rules for all index permutations with metric signature factors
- Simpler to compute CS terms from known structure than fully automate

### Attempted (Failed) Approaches
1. **Direct epsilon in Lagrangian:**
   ```mathematica
   (* Symbolic epsilon remains after decomposition *)
   CSTerm = epsilon[a,b,c] A[-a] CD[-b][A[-c]]
   (* Result: eps[{0,-cart3},{1,-cart3},{0,cart3}] ... *)
   (* Replacement rules complex and error-prone *)
   ```

2. **Field strength substitution:**
   ```mathematica
   (* Index mismatch - doesn't apply reliably *)
   F[-a,-b] -> CD[-a][A[-b]] - CD[-b][A[-a]]
   ```

## JSON Structure

### Manual JSON (Current)
File: `examples/data/chern_simons_3d.json`

```json
{
  "field": "A_0",
  "rhs": {
    "terms": [
      {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"},
      {"coefficient": 0.5, "operator": "gradient_x", "field": "A_2"},
      {"coefficient": -0.5, "operator": "gradient_y", "field": "A_1"}
    ]
  }
}
```

**Key points:**
- CS term decomposed into gradient operators (not a composite "curl")
- Cross-field references: A_0 couples to A_1 and A_2
- Python handles this naturally via existing gradient operators

### Why Manual JSON
ExportJSON.wl operator detection issues with 3D:
- Mixed Derivative[n,m] and Derivative[n,m,p] forms confuse parser
- Epsilon remnants in symbolic equations
- Cross-field detection works but coefficient extraction needs refinement

## Python Simulation

### Grid Setup (2D Spatial)
```python
grid = CartesianGrid(
    bounds=[(0, 50), (0, 50)],  # x and y spatial domain
    shape=[64, 64],
    periodic=True
)
```

### State Structure (6 Fields)
```python
# 3 vector components * (field + momentum)
state = FieldCollection([
    A_0, pi_0,  # Temporal component
    A_1, pi_1,  # x-spatial component
    A_2, pi_2   # y-spatial component
])
```

### Observed Physics
- Initial: Gaussian in A_1, others zero
- Evolution: Energy transfers to A_0 and A_2 via CS coupling
- Effect: Topological mass gives helical propagation patterns
- Verification: `max|A_0|` and `max|A_2|` grow from 0 as simulation runs

## Epsilon Tensor Component Values

For Minkowski signature (-,+,+) in 2+1D:

**Lower indices (covariant):**
```
ε_012 = -1  (includes √|det(g)| = 1 factor)
ε_021 = +1
ε_120 = -1
etc. (all permutations with sign)
```

**Upper indices (contravariant):**
```
ε^012 = +1  (raised with metric η^μν)
ε^021 = -1
etc.
```

**Critical:** Sign convention depends on metric signature. The built-in `epsiloneta3` from xAct uses standard conventions.

## Future Automation Plan

To fully automate CS term derivation:

1. **Define epsilon component rules:**
   ```mathematica
   epsilonRules = {
     epsiloneta3[{0,-cart3},{1,-cart3},{2,-cart3}] -> -1,
     (* ... all 6 permutations ... *)
   };
   ```

2. **Apply after ToBasis:**
   ```mathematica
   componentEq = ToBasis[chart][eom]
   componentEq = componentEq /. epsilonRules
   componentEq = Expand[componentEq]
   ```

3. **Enhance ExportJSON.wl:**
   - Detect curl patterns: `d_x A_y - d_y A_x`
   - Handle 3-arg Derivatives consistently
   - Map cross-derivative terms to gradient operators

4. **Test with full pipeline:**
   - Build CS Lagrangian with epsilon
   - VarD derives full EOM
   - Decompose and export to JSON
   - Verify coefficients and cross-field refs

## Related Examples

- **1+1D EM:** Maxwell without CS term (vector field baseline)
- **Coupled scalars:** Cross-field coupling pattern (scalar fields)
- **Chern-Simons:** Combines both: vector fields + cross-coupling + 2+1D

## References

- Script: `examples/chern_simons/chern_simons.wls`
- JSON: `examples/data/chern_simons_3d.json`
- Python: `examples/chern_simons/chern_simons_simulation.py`
- Plan: `/home/vscode/.claude/plans/valiant-orbiting-codd.md`
