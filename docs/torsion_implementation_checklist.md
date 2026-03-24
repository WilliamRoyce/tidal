# Torsion Implementation Checklist

## Stage 1: Pipeline Extensions

### 1a. Partial Antisymmetry Support
- [x] Modify `_field_definition()` to handle `symmetry = "antisymmetric_23"` (commit 419b758)
- [x] Add `symmetric_*` variant for completeness
- [x] Test: `massive_3form` example generates correct `Antisymmetric[{-a,-b,-c}]` (unchanged)
- [x] Test: partial `antisymmetric_23` generates correct `Antisymmetric[{-b,-c}]`
- [x] All 1531 tests pass
- [ ] Add dedicated unit test in `tests/test_cli.py`

### 1b. Torsion-Full CovD Support
- [x] Add `torsion = true` option to `[spacetime]` TOML config (commit ad98586)
- [x] Modify `_wls_spacetime()` to define CDT with `DefCovD[..., Torsion -> True]`
- [x] CDT name prefixing works for compound names (RicciScalarCDT, TorsionCDT, etc.)
- [x] CD prefixing NOT broken by CDT changes (negative lookahead regex)
- [x] All 1531 tests pass

### 1c. PGT Lagrangian Architecture (CRITICAL REVIEW FINDINGS)

**Investigated and resolved** (commit 0e81f84):

The original plan to use `ChangeCurvature + ChangeTorsion` was WRONG:
- `ChangeTorsion` converts TorsionCDT → Christoffel[CD,CDT] (contortion) — the WRONG direction
- After both operations, the Lagrangian contains only Christoffel[CD,CDT] — a derived quantity that DefTensorPerturbation CANNOT perturb
- xAct does not provide the inverse (contortion → torsion) natively

**Correct architecture**:
- User writes the Lagrangian in **already-decomposed form**: `RicciScalarCD[]` (Levi-Civita) + explicit `TorsionCDT[...]` terms
- Uses the identity R̃ = R^{LC} - 1/4 T² - 1/2 T²_bac + T²_trace + div (Shapiro 2002)
- `torsion = true` creates CDT and makes TorsionCDT available — it does NOT auto-decompose R̃
- DefTensorPerturbation on TorsionCDT correctly perturbs the torsion field
- The torsion tensor is NOT defined via `[[fields]]` — it comes from DefCovD

### 1d. Torsion Perturbation Registration
- [x] Implement automatic DefTensorPerturbation for TorsionCDT when torsion = true
- [x] Background torsion = 0 (flat Minkowski, no torsion)
- [x] xPert correctly perturbs TorsionCDT terms in the Lagrangian

### 1e. Run All Existing Tests
- [x] `uv run pytest tests/ -x -q` — all tests pass after each change
- [x] No regression in any rank-3 examples

## Stage 2: Graviton-Torsion Theory

### 2a. Documentation
- [x] Create `docs/tex/torsion.tex` with PGT framework, identity, architecture, references
- [x] Literature references (Shapiro 2002, Barker 2023, Hehl 1976, Hayashi 1979)

### 2b. Theory Config
- [x] Create `examples/graviton_torsion/theory.toml`
- [x] 4D Minkowski, torsion section, decomposed-form Lagrangian
- [x] Three T² invariants (α₁, α₂, α₃) + R̃² (b₅)
- [x] Linearization with automatic torsion perturbation registration
- [x] Plane-wave reduction (2+1D)

### 2c. Derivation
- [x] Component-level E-L derivation completes in ~5s (2+1D)
- [x] JSON: 15 fields (6 graviton + 9 torsion), 26 equations
- [x] Correct operator/coefficient structure verified
- [x] Hamiltonian: 14 terms with b₅ expressions, verified correct
- [x] Graviton-torsion mixing confirmed (derivative_3_x cross-coupling)

### 2d. Ostrogradsky Reduction
- [x] Automatic reduction of 4th-order-in-time equations on JSON load
- [x] 15 → 17 fields (2 auxiliary), all equations ≤ 2nd order

### 2e. Simulation & Validation
- [ ] **BLOCKED** by implicit d2_t/d3_t/mixed_T operators on 2nd-order fields (#165)
- [ ] Energy conservation
- [ ] Conversion measurement between h and T components
- [ ] Verify against analytical mixing matrix

## Stage 3: Component-Level E-L (General Pipeline — from torsion work)

- [x] ComponentEulerLagrange in ComponentDecompose.wl (commit 1b8e45d)
- [x] ReplaceHigherRankFieldComponents for antisymmetric tensors
- [x] Pure-time derivative filter in ExportJSON.wl
- [x] BuildMixedOperatorName for mixed time-space operators
- [x] derivative_3_x operator (3rd-order spatial, 5-point FD stencil)
- [x] Pipeline unification: component E-L default for ALL theories (commit 22e28a0)
- [x] 628 lines dead VarD code removed (commit 6406a4c)
- [x] CD shorthand optimisation extended to all theories (commit 21637ff)

## Key Findings Log

| Date | Finding | Impact |
|------|---------|--------|
| 2026-03-16 | ChangeTorsion goes wrong direction (T → K, not K → T) | Removed decomposition |
| 2026-03-16 | Christoffel[CD,CDT] is derived, not perturbable | Must keep TorsionCDT as fundamental |
| 2026-03-16 | TorsionCDT auto-created by DefCovD, not [[fields]] | Need special perturbation handling |
| 2026-03-20 | Abstract VarD + TraceBasisDummy: 45 indices → >77min | Motivated component E-L |
| 2026-03-21 | Component E-L: 5s vs 77min (900× speedup) | Now default for all theories |
| 2026-03-22 | R̃² produces 4th-order-in-time equations | Motivated Ostrogradsky reduction |
| 2026-03-22 | Ostrogradsky: 4th→2nd order automatic on JSON load | State: 15→17 fields |
| 2026-03-22 | Implicit d2_t cross-references block simulation | Issue #165, solver-side fix needed |
