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
- Group by subsystem: Wolfram/xAct, Python/py-pde, Pipeline Integration
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

## Common Python/py-pde Issues

### py-pde Type Annotation Patterns

**Symptoms:** pyright/mypy errors on py-pde method returns (`gradient()`, `laplace()`, arithmetic operators)

**Cause:** py-pde's type stubs don't fully describe runtime return types. Field arithmetic and gradient operations return generic base types, not `ScalarField`.

**Solutions:**

1. Use `TYPE_CHECKING` imports for py-pde types to avoid runtime overhead:
   ```python
   from typing import TYPE_CHECKING, cast
   if TYPE_CHECKING:
       from pde import CartesianGrid, FieldCollection, ScalarField
   ```

2. Cast after arithmetic and gradient operations:
   ```python
   result = phi - psi
   result = cast("ScalarField", result)  # py-pde returns DataFieldBase at type level

   grad = phi.gradient(bc=cast("Any", bc))
   grad = cast("FieldCollection", grad)  # gradient returns VectorField at runtime
   ```

3. Always add runtime `isinstance` checks before casts — fail-fast on type mismatches:
   ```python
   if not isinstance(dpi_dt, ScalarField):
       msg = "dpi_dt computed non-ScalarField result"
       raise TypeError(msg)
   ```

**Don't:** Suppress type errors with `# type: ignore` — use casts so the type checker still validates downstream usage.

### Boundary Condition Gotcha: Gradient Chaining

**Symptoms:** `gradient().gradient()` (for directional second derivatives) gives wrong results or silent numerical drift on periodic grids

**Cause:** py-pde requires explicit `"auto_periodic_neumann"` boundary condition string for gradient chaining. Passing `None` for periodic grids silently breaks the computation.

**Solutions:**

1. For periodic grids, always use:
   ```python
   bc = "auto_periodic_neumann"
   d2_phi_dx2 = phi.gradient(bc)[0].gradient(bc)[0]  # correct ∂²φ/∂x²
   ```

2. For non-periodic grids, use `"derivative"` (Neumann, zero flux):
   ```python
   bc = "derivative"
   ```

3. For mixed periodic/non-periodic (per-axis), check `grid.periodic`:
   ```python
   periodic = getattr(grid, "periodic", None)
   if isinstance(periodic, Sequence):
       bc = "auto_periodic_neumann" if any(periodic) else "derivative"
   ```

See `infer_bc_from_grid()` in `tidal/utils.py` which encapsulates this logic.

### Grid Coordinate Version Compatibility

**Symptoms:** `IndexError` or wrong-shaped arrays when accessing `grid.cell_coords`

**Cause:** Different py-pde versions return cell coordinates in different formats:
- Grid-shaped: `(*grid.shape, dim)` — newer versions
- Flattened: `(N_cells, dim)` — older versions

**Solution:** Always normalize before use:
```python
coords = cast("np.ndarray", grid.cell_coords)
if coords.ndim == grid.dim + 1:  # grid-shaped
    coords = coords.reshape(-1, grid.dim)  # flatten
# Now coords is always (N_cells, dim)
```

### Numba Backend Fallback

**Symptoms:** `NotImplementedError` when solving with `backend="numba"`

**Cause:** Not all PDE classes implement `_make_pde_rhs_numba`. Custom PDEs using `evolution_rate()` override instead of expression-based definitions can't auto-compile to Numba.

**Solutions:**
1. Check before attempting: `hasattr(pde, "_make_pde_rhs_numba")`
2. Fallback pattern:
   ```python
   try:
       result = pde.solve(state, ..., backend="numba")
   except NotImplementedError:
       result = pde.solve(state, ..., backend="numpy")
   ```

### Spatial Coefficient Freezing for Numba JIT

**Symptoms:** Numba JIT functions see stale or corrupted spatial coefficient data

**Cause:** Numba closures capture references, not copies. If the original array is modified, the JIT function sees the new data.

**Solution:** Copy spatial data before creating the JIT closure:
```python
m2_data = self.m2_field.data.copy()  # freeze before JIT
laplace = state.grid.make_operator("laplace", bc)  # create operator outside JIT

@jit
def pde_rhs(state_data, t):
    return laplace(state_data) - m2_data * state_data  # uses frozen copy
```

### Solver Return Type Polymorphism

**Symptoms:** `TypeError` or `AttributeError` accessing solve result

**Cause:** `pde.solve()` can return either `FieldCollection` directly or `tuple[FieldCollection | None, dict]`.

**Solution:** Use `normalize_solve_result()` from `tidal/utils.py`:
```python
from tidal.utils import normalize_solve_result
raw = pde.solve(state, t_range=t_end, ...)
result = normalize_solve_result(raw)  # always FieldCollection, raises if None
```

### "Unknown operator: X"

**Symptoms:** `ValueError` when building PDE

**Cause:** JSON references operator not implemented in `pde_builder.py`

**Available operators:**
- `identity`, `laplacian`
- `gradient_x`, `gradient_y`, `gradient_z`

**Fix:** Add operator to `_get_operator` method in `pde_builder.py`

### Grid Dimension Mismatch

**Symptoms:** Error about field dimensions vs grid dimensions

**Cause:** Using 1D grid for 2D spatial problem (or vice versa)

**Check:**
- 1+1D: `CartesianGrid(bounds=[(0,100)], shape=[256])`  # 1 spatial dimension
- 2+1D: `CartesianGrid(bounds=[(0,50),(0,50)], shape=[64,64])`  # 2 spatial dimensions

### State Size Mismatch

**Symptoms:** `AssertionError` about state size

**Cause:** State has wrong number of fields

**Correct sizes:**
- N components: 2N fields (N fields + N momenta)
- Example: 3 vector components → 6 total fields

**Fix:**
```python
# For 3 components
state = FieldCollection([A_0, pi_0, A_1, pi_1, A_2, pi_2])
assert len(state) == 6
```

### Storage Tracker Error

**Symptoms:** `TypeError: unexpected keyword argument 'interval'`

**Cause:** py-pde API uses positional argument, not keyword

**Fix:**
```python
# BAD
tracker=storage.tracker(interval=0.5)

# GOOD
tracker=storage.tracker(0.5)
```

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
- [ ] Run all examples: `wolframscript -file examples/*/*.wls`
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
