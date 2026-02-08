# Project Memory: Torsion-Gertsenshtein Lagrangian-to-PDE Pipeline

## 📝 Document Maintenance

**IMPORTANT:** This file should be updated throughout the project lifecycle to capture learnings and patterns.

**When to update this file:**
- After implementing a new physics example (scalar, vector, tensor fields)
- When discovering new xAct/Wolfram patterns or pitfalls
- After generalizing infrastructure (dimension support, operator handling)
- When solving tricky bugs that reveal architectural insights
- After making critical design decisions that affect future work

**What to update:**
- Add new patterns to relevant sections (xAct, Pipeline, Python)
- Document workarounds and their reasons in "Known Issues & Solutions"
- Update "Example Implementations" with new examples and their key features
- Expand "Quick Reference" with new file locations and key functions
- Add new sections for major features (e.g., gauge fixing, nonlinear terms)

**Companion files:**
- `troubleshooting.md` - Error patterns and debugging techniques
- `chern-simons-notes.md` - Example-specific implementation details
- Create similar notes files for complex examples (e.g., `yang-mills-notes.md`)

---

## Overview
This project implements a symbolic physics pipeline: Lagrangian (xAct/Mathematica) → JSON → PDE simulation (py-pde). Users should NEVER manually hardcode equations - all PDEs must derive from Lagrangian via symbolic computation.

## Critical Architecture Decisions

### Multi-Field Coupling (WORKING)
- **Cross-field terms:** Use gradient operators in JSON (e.g., `gradient_x(A_2)` in A_0 equation)
- **Field transformation:** `DecomposeToComponents` accepts `additionalFields` parameter to transform all coupled fields to coordinate form
- **JSON detection:** `IdentifyMultiFieldTerm` uses case-insensitive symbol matching and function head extraction

### Dimension Handling (GENERALIZED)
- **Dynamic dimension:** `ComponentDecompose.wl` now uses `dim = Length[ScalarsOfChart[chart]]` instead of hardcoded `dim = 2`
- **Derivative conversion:** `CommonUtilities.wl` handles both 2-arg (1+1D) and 3-arg (2+1D) Derivative forms
- **Supports:** 1+1D (t,x), 2+1D (t,x,y), extensible to 3+1D

## xAct/Wolfram Patterns

### Symbol Management
- **Always check existence:** `If[!xTensorQ[M2], DefManifold[M2, ...]]`
- **Use shared symbols:** M2, M3, eta, CD for examples to avoid conflicts
- **Avoid random names:** Kernel caching makes RandomInteger-based names fail

### Levi-Civita / Epsilon Tensors
- **Built-in epsilon:** `DefMetric` automatically creates `epsiloneta3[-a,-b,-c]`
- **Manual components:** Epsilon tensor values need explicit evaluation in coordinates
- **Current limitation:** Full symbolic epsilon derivation not automated - use manual addition for CS-like terms
- **Workaround pattern:** Derive Maxwell part symbolically, add topological terms manually with known structure

### Field Strength Tensors
- **Expand before decompose:** `csF[-a,-b]` must be replaced with `CD[-a][csA[-b]] - CD[-b][csA[-a]]` BEFORE `DecomposeToComponents`
- **Index matching:** Abstract index patterns in rules don't auto-apply - construct Lagrangian directly in terms of derivatives of A

### Lagrangian Construction
```mathematica
(* GOOD: Direct construction *)
L = -1/2 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]]

(* AVOID: Field strength with substitution *)
L = -1/4 F[-a,-b] F[a,b]  (* Then F -> CD[A] - CD[A] doesn't reliably apply *)
```

## Pipeline Module Integration

### Loading Pattern
```mathematica
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];
```

### Cross-Field Decomposition
```mathematica
(* For coupled fields phi and chi *)
phiComponents = DecomposeToComponents[eomPhi, phi[], cart, {chi[]}];
chiComponents = DecomposeToComponents[eomChi, chi[], cart, {phi[]}];
```

### JSON Export
```mathematica
(* Multi-field format *)
fieldEquations = {{"phi_0", eqPhi}, {"chi_0", eqChi}};
jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];
```

## Known Issues & Solutions

### Issue: Component equations return 0
**Cause:** Field strength tensor not expanded before decomposition
**Solution:** Replace `csF` with derivative expression before calling `DecomposeToComponents`

### Issue: Cross-field terms not detected
**Cause:** Other fields not transformed to coordinate form (e.g., `cplChi[]` vs `cplChi0[t,x]`)
**Solution:** Pass all coupled fields to `DecomposeToComponents(..., additionalFields={...})`

### Issue: Coefficients all show as 1.0
**Cause:** Pattern matching in `IdentifyMultiFieldTerm` not finding field symbols
**Solution:** Extract function heads (not just symbols), use case-insensitive matching

### Issue: "Symbol X already used as manifold"
**Cause:** xAct kernel caching across script runs
**Solution:** Use `xTensorQ[X]` check before `DefManifold`, reuse standard symbols (M2, M3)

### Issue: Mixed 2-arg and 3-arg Derivatives in 2+1D
**Cause:** `ConvertCDToDerivatives` creates inconsistent forms
**Status:** Known limitation - JSON parser handles but can confuse operator detection
**Workaround:** Manual JSON creation for complex 2+1D examples (see chern_simons_3d.json)

## Example Implementations

### Coupled Scalars (WORKING)
- **File:** `examples/coupled_scalars/coupled_scalars.wls`
- **Features:** Cross-field coupling, mass matrix, energy transfer
- **JSON:** Auto-generated via `BuildMultiFieldJSONStructure`
- **Key:** Pass `{chi[]}` when decomposing phi, `{phi[]}` when decomposing chi

### Chern-Simons 2+1D (WORKING - AUTOMATED)
- **File:** `examples/chern_simons/chern_simons.wls`
- **Features:** 3D manifold, full Maxwell-CS symbolic derivation, auto epsilon evaluation
- **JSON:** Auto-generated with gradient_x, gradient_y operators
- **Key:** Use `epsiloneta3[a, b, c]` in Lagrangian, pipeline handles rest
- **New functions:** `EvaluateEpsilonComponents`, `IdentifyGradientDirection`
- **Limitation:** Mixed time-space derivatives in 2+1D need better handling

## Testing Guidelines

### Regression Testing
- Run existing 1+1D examples after Wolfram changes: `wolframscript -file examples/*/*.wls`
- Run pytest suite: `uv run pytest tests/` (expect 496 Python tests + ~100 Wolfram tests passing)
- Check for dimension hardcoding in error messages

### Verification Pattern
1. Wolfram derivation produces N component equations
2. JSON has correct dimension, signature, coordinates
3. JSON terms show proper cross-field references and coefficients
4. Python simulation runs without errors
5. Physical behavior matches expectations (e.g., energy transfer)

## Python Side (py-pde)

### Operators Supported
- `identity`: field itself
- `laplacian`: spatial Laplacian (works in 1D, 2D, 3D)
- `gradient_x`, `gradient_y`, `gradient_z`: directional gradients
- All operators support cross-field application

### Grid Requirements
- 1+1D: CartesianGrid with 1D bounds (t is time)
- 2+1D: CartesianGrid with 2D bounds (t is time, x and y are spatial)
- Always use periodic BCs for wave equations

## Completed Improvements

1. ~~**Automate epsilon tensor handling:**~~ ✅ DONE - `EvaluateEpsilonComponents` in CommonUtilities.wl
2. ~~**Generalize to 3+1D:**~~ ✅ DONE - Full 4D spacetime with 3 spatial dimensions (gravitational_waves, scalar_field_3d)
3. ~~**Higher-rank tensors:**~~ ✅ DONE - Phase 13: Rank 3+ tensor support (Issue #70)
4. ~~**Time derivative detection in 2+1D:**~~ ✅ DONE - Issue #79: Mixed time-spatial derivatives
5. ~~**Unified derivative classification:**~~ ✅ DONE - Issue #85: Dimension-agnostic `ExtractDerivativeProfile`

## Future Improvements Needed

1. **Gauge fixing:** Automate Lorenz/Coulomb gauge application in Wolfram
2. **Clean JSON output for 2+1D:** Remove spurious cross-Laplacian terms from mixed derivatives (low priority)
3. **Nonlinear extensions:** Beyond linear perturbation theory
4. **Continuous Integration:** GitHub Actions for automated Wolfram test execution

## Quick Reference

### File Locations
- Wolfram modules: `torsion_gertsenshtein/wolfram/`
- Python pipeline: `torsion_gertsenshtein/symbolic/`
- Examples: `examples/{scalar_field,electromagnetic,coupled_scalars,chern_simons}/`
- Generated JSON: `examples/data/*.json`

### Key Functions
- `DecomposeToComponents[eom, field, chart, additionalFields]` - Component extraction
- `BuildMultiFieldJSONStructure[fieldEquations, metadata]` - JSON generation
- `build_pde_from_json(json_path)` - Python PDE construction
- `IdentifyMultiFieldTerm[term, currentField, allFields, metadata]` - Operator detection

See detailed notes in: `/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/`
