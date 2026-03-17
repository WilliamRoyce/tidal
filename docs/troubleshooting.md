# Troubleshooting Guide

## 📝 Maintaining This Guide

**CRITICAL:** Update this file immediately when you encounter and solve new errors. Future work depends on these patterns being documented.

**When to add entries:**

- **Always** when you hit a new error that takes more than 10 minutes to debug
- When existing error patterns manifest in new ways
- After discovering a non-obvious cause for a common symptom
- When you find a better solution to an existing problem
- After adding new features that create new failure modes

**How to structure entries:**

```
### Error Title (Brief, Searchable)

**Symptoms:** What the user sees (error messages, wrong output)
**Cause:** Root cause explanation
**Solutions:** Numbered steps to fix
**Don't:** Anti-patterns that seem like they'd work but don't
```

**Organization:**

- Group by subsystem: Wolfram/xAct, Python solver, Pipeline Integration
- Keep "Debugging Techniques" and "Verification Checklist" sections at the end
- Cross-reference with `MEMORY.md` for architectural context
- Link to example-specific notes (like `chern-simons-notes.md`) for complex cases

**Pruning:** Remove entries if:

- The underlying code changed and the error can't happen anymore
- Better solutions made the workaround obsolete (note the improvement in MEMORY.md)

---

## Common Wolfram/xAct Issues

### "Symbol X is already used as a manifold"

**Symptoms:** Error when running script multiple times in same kernel session

**Cause:** xAct caches tensor definitions in kernel memory

**Solutions:**

1. Always check before defining:

   ```mathematica
   If[!xTensorQ[M2], DefManifold[M2, 3, {a,b,c,d,e,f}]]
   ```

2. Use standard shared symbols across examples:
   - Manifolds: M2 (1+1D), M3 (2+1D), M4 (3+1D)
   - Metrics: eta, eta3, eta4
   - Covariant derivatives: CD, CD3, CD4
   - Charts: cart, cart3, cart4

3. Restart kernel if needed: `Quit[]` or `wolframscript -file` (fresh process each time)

**Don't:** Use `RandomInteger` or timestamp-based symbol names - kernel caching makes them fail

### Component Equations Return All Zeros

**Symptoms:** After `DecomposeToComponents`, all equations show `0`

**Likely causes:**

1. **Field strength not expanded:**

   ```mathematica
   (* BAD *)
   eom = VarD[A[-a], CD][L]  (* L contains F[-a,-b] *)
   components = DecomposeToComponents[eom, A[-a], cart]  (* F not expanded *)

   (* GOOD *)
   eom = VarD[A[-a], CD][L /. F[-a,-b] -> CD[-a][A[-b]] - CD[-b][A[-a]]]
   (* OR construct L directly in terms of CD[A] *)
   ```

2. **Missing field in ToBasis conversion:**
   - Check that field symbols appear in componentEq after `ToBasis[chart][eom]`
   - If abstract field not converted, check field definition and chart compatibility

### Cross-Field Terms Not Detected

**Symptoms:** JSON shows all terms referencing same field, no cross-coupling

**Cause:** Other fields not transformed to coordinate form

**Example problem:**

```mathematica
(* Chi appears as cplChi[] not cplChi0[t,x] *)
phiEq = -0.5*cplChi[] - 1.0*cplPhi0[t,x] + Derivative[...][cplPhi0][t,x]
```

**Solution:**

```mathematica
(* Pass all coupled fields *)
phiComponents = DecomposeToComponents[eomPhi, phi[], cart, {chi[]}]
chiComponents = DecomposeToComponents[eomChi, chi[], cart, {phi[]}]
```

**Verification:**

- After decomposition, print equations
- All field symbols should have coordinate arguments: `field0[t,x]` or `field0[t,x,y]`
- No bare field symbols like `field[]`

### JSON Coefficients All Show 1.0

**Symptoms:** Correct operator/field but coefficient extraction fails

**Cause:** Pattern matching in `ExtractNumericCoefficient` not finding field symbols

**Debug steps:**

1. Check `IdentifyMultiFieldTerm` function head extraction
2. Verify field names match: "phi_0" → "phi" → "Phi" (case variations)
3. Print intermediate term structure to see actual symbols

**Recent fix:** Use `ToLowerCase` for field base name matching (case-insensitive)

### Multiline Lagrangian Parsed Incorrectly

**Symptoms:** Second field's terms evaluate to zero

**Cause:** Mathematica multiline without explicit `+`

**Solution:**

```mathematica
(* BAD *)
L = term1
    term2  (* Treated as separate expression *)

(* GOOD *)
L = term1 +
    term2  (* Explicit continuation *)

(* OR *)
L = (
  term1 +
  term2
)
```

**Example from coupled scalars:** Required `+` before `(-1/2 CD[-a][chi[]]...)`

### Epsilon Tensor Not Evaluating to Numbers

**Symptoms:** Epsilon tensor like `epsiloneta3[{0,-cart3},{1,-cart3},{2,cart3}]` remains in output

**Cause:** Pattern mismatch in `EvaluateEpsilonComponents`

**Common issues:**

1. Wrong chart name (must match exactly)
2. Index sign pattern not handled (mixed up/down indices)
3. Function not being called in pipeline

**Solution:**

```mathematica
(* Verify epsilon is being evaluated *)
testExpr = epsiloneta3[{0, -cart3}, {1, -cart3}, {2, cart3}];
result = EvaluateEpsilonComponents[testExpr, cart3];
Print["Result: ", result];  (* Should be a number, not symbolic *)
```

**Debug:** If epsilon remains symbolic, check:

- The chart variable matches exactly (e.g., `cart3` not `cart`)
- The epsilon tensor name matches pattern (must contain "epsilon")
- Both covariant (`-chart`) and contravariant (`chart`) cases are handled

### Mixed Time-Space Derivatives in 2+1D

**Symptoms:** Extra "laplacian" terms on cross-fields in JSON

**Cause:** `Derivative[1, 1, 0][f][t,x,y]` (d_t d_x f) classified as second-order

**Current status:** Known limitation - the time derivative detection pattern is 1+1D only

**Workaround:** For now, the extra terms don't break simulation but may cause numerical artifacts

## Common Python Solver Issues

> **Note:** TIDAL's solver layer uses SUNDIALS IDA/CVODE + native numpy operators (see `tidal/solver/`). The former py-pde-based architecture was replaced in February 2026.

### "Unknown operator: X"

**Symptoms:** `ValueError` when loading or running a simulation

**Cause:** JSON spec references an operator not in `OPERATOR_REGISTRY` (`tidal/solver/operators.py`)

**Available operators:** `identity`, `laplacian`, `laplacian_x`, `laplacian_y`, `laplacian_z`, `gradient_x`, `gradient_y`, `gradient_z`, `cross_derivative_xy`, `cross_derivative_xz`, `cross_derivative_yz`, `biharmonic`, `first_derivative_t`, `mixed_T_S1_S2_...`

**Fix:** Check `OPERATOR_REGISTRY` in `tidal/solver/operators.py`. If the operator is valid but missing, add it there.

### Grid Dimension Mismatch

**Symptoms:** Error about field dimensions vs grid dimensions

**Cause:** Using `--grid-shape` with wrong number of dimensions for the JSON spec's spatial dimension

**Check:** Match grid shape to the JSON spec's `spacetime.dimension - 1` spatial dimensions:

- 1+1D specs: `--grid-shape 256` (1 spatial dimension)
- 2+1D specs: `--grid-shape 64,64` (2 spatial dimensions)
- 3+1D specs: `--grid-shape 32,32,32` (3 spatial dimensions)

### State Size Mismatch

**Symptoms:** Error about unexpected state vector size

**Cause:** State has wrong number of slots for the equation system

**How slot counts work** (via `StateLayout` in `tidal/solver/state.py`):

- 2nd-order fields: 2 slots each (field + velocity, e.g., `phi_0` + `v_phi_0`)
- 1st-order fields: 1 slot each (field only)
- Constraint fields: 1 slot each (solved algebraically)
- Example: 3 vector components (2nd order) → 6 total slots

### Non-Periodic Coefficient Warning / Error

**Symptoms:** Warning or `ValueError` at simulation start:
`"Position-dependent coefficient '...' has N% jump at the periodic boundary along x (left=..., right=..., leak_metric=...)"`

**Cause:** Position-dependent coefficients (e.g., background magnetic field `B(x)`) must be
continuous across periodic boundaries. If they aren't, the integration-by-parts identity
fails and causes O(1) energy non-conservation that does **not** improve with finer grids.

The check uses a *leak metric* = `(rel_jump) * (boundary_significance)` to estimate the
actual energy leak magnitude. Warnings fire when the metric exceeds 0.01; errors when it
exceeds 0.25.

**Solutions:**

1. **Large domain**: Use a domain large enough that the non-periodic coefficient (e.g.,
   dipolar `B ~ 1/r^3`) is negligible at both boundaries. Terminate the simulation before
   waves reach the edges.
2. **Localised profile**: Use a windowed/Gaussian profile that naturally goes to the same
   value at both boundaries (e.g., `B(x) = Bpeak * exp(-(x-x0)^2/R^2)`).
3. **Non-periodic BCs**: Switch to Dirichlet or Neumann BCs (`--bc neumann`).
4. **False positive?**: If the warning fires but energy conservation is actually fine
   (check with `tidal measure --what conservation`), the coefficient boundary values are
   small enough that the leak is negligible. The leak metric should be very small in this
   case — if it isn't, the warning is genuine.

### Constraint Solver Failures

**Symptoms:** IDA fails to converge, or constraint fields have NaN values

**Cause:** Algebraic constraint equations may have singular operators (e.g., Laplacian with periodic BCs has a null space for constant functions)

**Solutions:**

1. The three-tier constraint pre-solve (`tidal/solver/constraint_solve.py`) handles most cases automatically
2. For periodic BCs with pure-Laplacian constraints, gauge regularisation pins the zero mode
3. Check `--ic` settings — constraint ICs must be consistent with dynamical field ICs

## Debugging Techniques

### Wolfram Side

1. **Print intermediate steps:**

   ```mathematica
   Print["After VarD: ", eom];
   Print["After ToBasis: ", ToBasis[chart][eom]];
   Print["After metric eval: ", EvaluateMinkowskiMetric[expr, chart]];
   ```

2. **Check field transformation:**

   ```mathematica
   (* Should see coordSyms = {t[], x[]} or {t[], x[], y[]} *)
   coordSyms = GetCoordinateSymbols[chart]
   Print[coordSyms]
   ```

3. **Verify dimension:**

   ```mathematica
   dim = Length[ScalarsOfChart[chart]]
   Print["Dimension: ", dim]
   ```

4. **Test ToCanonical:**
   ```mathematica
   (* Simplify before and after to see if indices contract properly *)
   Print["Before: ", expr];
   expr = ToCanonical[expr];
   expr = ContractMetric[expr];
   Print["After: ", expr];
   ```

### Python Side

1. **Validate JSON load:**

   ```python
   spec = load_equation_system(json_path)
   print(f"Dimension: {spec.dimension}")
   print(f"Components: {spec.component_names}")
   for eq in spec.equations:
       print(f"{eq.field_name}: {len(eq.rhs_terms)} terms")
   ```

2. **Check state structure:**

   ```python
   print(f"State fields: {len(state)}")
   for i, field in enumerate(state):
       print(f"  {i}: {field.label}, shape={field.data.shape}")
   ```

3. **Monitor evolution:**
   ```python
   # Add print inside evolution_rate or use callback tracker
   def check_rates(state):
       print(f"Max field value: {max(np.max(np.abs(f.data)) for f in state)}")
   ```

## Verification Checklist

### After Wolfram Changes

- [ ] Run all examples: `./scripts/run_examples.sh` (or `tidal derive examples/*/theory.toml`)
- [ ] Check JSON dimension matches spacetime (2 for 1+1D, 3 for 2+1D)
- [ ] Verify component count matches field rank × dimensions
- [ ] Spot-check coefficients in JSON (shouldn't all be 1.0)
- [ ] Confirm cross-field terms if applicable

### After Python Changes

- [ ] Run pytest: `uv run pytest tests/`
- [ ] Test JSON loading for all examples
- [ ] Run at least one simulation end-to-end
- [ ] Check output plots if visualization enabled
- [ ] Verify energy conservation (if applicable to physics)

### After Pipeline Changes

- [ ] Regression test: coupled_scalars still works
- [ ] Forward test: new example runs
- [ ] Cross-test: modify old example to use new feature
- [ ] Documentation: update MEMORY.md with new patterns
