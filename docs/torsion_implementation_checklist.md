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

**Remaining question**: How to register TorsionCDT for perturbation. Options:
- (a) Auto-register when `torsion = true` + linearization present (cleanest)
- (b) Special `[[linearization.torsion_perturbation]]` TOML section
- (c) Reference via `[[linearization.matter_perturbations]]` with field = "TorsionCDT" (needs special handling since it's not in [[fields]])

### 1d. Torsion Perturbation Registration
- [ ] Implement automatic DefTensorPerturbation for TorsionCDT when torsion = true
- [ ] Background torsion = 0 (flat Minkowski, no torsion)
- [ ] Test: dry-run generates correct DefTensorPerturbation code
- [ ] Test: xPert correctly perturbs TorsionCDT terms in the Lagrangian

### 1e. Run All Existing Tests
- [x] `uv run pytest tests/ -x -q` — all 1531 tests pass after each change
- [ ] No regression in any rank-3 examples

## Stage 2: Graviton-Torsion Theory

### 2a. Documentation
- [ ] Create `docs/torsion.md` with:
  - PGT framework (tetrad, connection, torsion as field strength)
  - The R̃ = R^{LC} + T² identity (Shapiro 2002, eq 2.17)
  - Why ChangeCurvature/ChangeTorsion doesn't work (and what does)
  - How to write the TOML Lagrangian
  - xAct torsion functions and their directions
  - Literature references

### 2b. Theory Config
- [ ] Create `examples/graviton_torsion/theory.toml`:
  - 4D Minkowski, `torsion = true`
  - Lagrangian: `(1/κ²) RicciScalarCD[] + α₁ TorsionCDT² + α₂ TorsionCDT²_bac + α₃ TorsionCDT²_trace`
  - Linearization: metric perturbation h + torsion perturbation (auto-registered)
  - Plane-wave reduction
  - TT gauge for graviton only

### 2c. Derivation
- [ ] `tidal derive` completes successfully
- [ ] Inspect JSON: fields, operators, coupling structure
- [ ] Check Einstein-Cartan limit (α_I = 0)

### 2d. Simulation & Validation
- [ ] Energy conservation
- [ ] Conversion measurement between h and T components
- [ ] Verify against analytical mixing matrix

## Key Findings Log

| Date | Finding | Impact |
|------|---------|--------|
| 2026-03-16 | ChangeTorsion goes wrong direction (T → K, not K → T) | Removed decomposition |
| 2026-03-16 | Christoffel[CD,CDT] is derived, not perturbable | Must keep TorsionCDT as fundamental |
| 2026-03-16 | TorsionCDT auto-created by DefCovD, not [[fields]] | Need special perturbation handling |
