# Torsion Implementation Checklist

## Stage 1: Pipeline Extensions

### 1a. Partial Antisymmetry Support
- [ ] Modify `_field_definition()` in `tidal/cli/_derive.py` to handle `symmetry = "antisymmetric_23"` etc.
- [ ] Add `symmetric_*` variant for completeness
- [ ] Test: existing `massive_3form` example (fully antisymmetric) still derives correctly
- [ ] Test: new rank-3 partial-antisymmetric field generates correct xAct `Antisymmetric[{-b,-c}]`
- [ ] Test: dry-run derive with partial antisymmetry produces valid Wolfram code
- [ ] Add unit test in `tests/test_cli.py` for partial antisymmetry code generation

### 1b. Torsion-Full CovD Support
- [ ] Add `torsion = true` option to `[spacetime]` TOML config
- [ ] Modify `_wls_spacetime()` to define CDT with `DefCovD[..., Torsion -> True]`
- [ ] Generate `ChangeCurvature`/`ChangeTorsion` decomposition code before xPert
- [ ] The Lagrangian can reference `RicciScalarCDT[]` and `TorsionCDT[...]`
- [ ] Test: verify `RicciScalarCDT[]` evaluates correctly in Wolfram
- [ ] Test: verify decomposition R̃ → R^{LC} + T² terms gives correct result
- [ ] Add unit test for torsion CovD TOML parsing

### 1c. Torsion as Matter Perturbation
- [ ] Verify `[[linearization.matter_perturbations]]` works for rank-3 tensor with partial antisymmetry
- [ ] The torsion field T background = 0 (flat Minkowski, no torsion)
- [ ] `DefTensorPerturbation` creates perturbation of `TorsionCDT`
- [ ] Test: dry-run derive with torsion perturbation

### 1d. Run All Existing Tests
- [ ] `uv run pytest tests/ -x -q` — all 1531+ tests pass
- [ ] No regression in massive_3form or any rank-3 example

## Stage 2: Graviton-Torsion Theory

### 2a. Documentation
- [ ] Create `docs/torsion.md` with:
  - PGT framework (tetrad, connection, torsion as field strength)
  - The R̃ = R^{LC} + T² identity
  - How the TOML Lagrangian is constructed
  - Torsion irreducible decomposition (tensor/trace/axial)
  - xAct torsion functions used
  - Literature references

### 2b. Theory Config
- [ ] Create `examples/graviton_torsion/theory.toml` with:
  - 4D Minkowski spacetime with `torsion = true`
  - Fields: h (rank-2 symmetric), T (rank-3, antisymmetric_23) via torsion CovD
  - Lagrangian: (1/κ²)R̃ + α₁T²_{abc} + α₂T²_{bac} + α₃T²_trace
  - Constants: kappa, alpha1, alpha2, alpha3
  - Linearization with metric perturbation + torsion perturbation
  - Plane-wave reduction
  - TT gauge for graviton (NOT for torsion)

### 2c. Derivation
- [ ] `tidal derive examples/graviton_torsion/theory.toml` completes successfully
- [ ] Inspect JSON: correct fields, operators, coupling structure
- [ ] Identify which torsion components survive constraint elimination
- [ ] Check Einstein-Cartan limit (α_I = 0): all torsion components eliminated

### 2d. Simulation & Validation
- [ ] Simulate with modal solver or CVODE
- [ ] Energy conservation |dE/E| < 1e-6
- [ ] Conversion measurement between graviton and torsion modes
- [ ] Verify at α_I = 0: no conversion (Einstein-Cartan limit, torsion = 0)
- [ ] Verify at α_I ≠ 0: Rabi-like oscillation between h and T components

## Stage 3: Sweeps & Benchmarks (future)
- [ ] Sweep α_I parameters
- [ ] Compare with analytical mixing matrix eigenvalues
- [ ] Create sweep scripts with analytical overlays
- [ ] Create `examples/graviton_torsion/run.sh`

## Stage 4: Add EM Sector (future)
- [ ] Add Maxwell Lagrangian to theory
- [ ] Non-minimal torsion-EM coupling
- [ ] 3-body mixing: graviton-torsion-photon
