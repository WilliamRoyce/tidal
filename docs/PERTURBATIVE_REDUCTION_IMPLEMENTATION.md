# Perturbative Reduction Implementation Tracking (v6)

Tracks the v6 iterative-numerical perturbative expansion. Supersedes the v5 JLM symbolic reduction approach (deleted). Full plan: `/home/vscode/.claude/plans/flickering-gathering-orbit.md`. GitHub issue: [#267](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/267).

**Architecture (v6, 2026-04-19)**: Treat small-parameter corrections (b₅R̃², ρ·EH, etc.) as sources for a second solver pass atop the unmodified base equations. Theory-agnostic, ghost-free, no classification needed. Phase C in Wolfram reduces to order-tagging in JSON; all perturbative work happens at simulate time in Python.

**Literature** (both in `literature/`): `gr-qc_9211002/` (Parker-Simon 1993), `2505.00082/` (FKY 2025).

## Status Legend

- [ ] Not started
- [~] In progress
- [x] Completed
- [!] Blocked / needs attention

---

## Stage 0: Extend Post-E-L Transforms for Dual Application — COMPLETE

**Kept from v5**: Independent cleanup. Gauge fixing, plane-wave reduction, TT substitution extended to modify `lagComp` as well as `fieldEquations`. See prior tracking log for details.

### Tasks

- [x] **0.1** Extend `_wls_gauge_fixing_type_b(ctx)` in `tidal/cli/_derive.py` to also modify `lagComp`
- [x] **0.2** Extend `_wls_plane_wave_reduction_equations(ctx)` to `lagComp`
- [x] **0.3** Extend `_wls_plane_wave_field_elimination(ctx)` to `lagComp`
- [x] **0.4** Extend `_wls_plane_wave_coordinate_evaluation(ctx)` to `lagComp`
- [x] **0.4b** Extend `_tt_traceless_substitution` to `lagComp`
- [x] **0.5** Smoke-test on `gertsenshtein` theory (WLS parses, CLI tests pass)

---

## Stage 1: Order-tagging in ExportJSON.wl + delete PerturbativeReduction.wl — COMPLETE

**Goal**: Annotate every OperatorTerm with `order_in_eps` = total exponent of small parameters in its symbolic coefficient. Delete the v5 JLM module.

**Estimated**: 3 days. **Actual**: 1 day (2026-04-19).

### Tasks

- [x] **1.1** Extended `BuildTermResult` in `tidal/wolfram/ExportJSON.wl:862` with `orderInEps_:0` parameter; emits `"order_in_eps"` field when > 0.
- [x] **1.2** Added `ComputeOrderInEps[expr, paramsList]` helper in ExportJSON.wl:886. Uses `Max[Total[Exponent[term, p]] for p in params]` across Expand'd monomials. Clamps negative exponents to 0.
- [x] **1.3** Threaded `smallParams` through `ExtractTermCoefficient` → `IdentifyMultiFieldTerm` → `BuildTermResult`. `ExtractTermCoefficient` now returns a 5-tuple including `orderInEps`.
- [x] **1.4** `_derive.py` injects `"small_parameters" -> {param_sym_list}` and `"perturbation_order" -> N` into Wolfram metadata Association when `ctx.perturbative_reduction` is set. `EquationToJSONMultiField` reads these via `Lookup[metadata, "small_parameters", {}]`. `BuildMultiFieldJSONStructure` also emits a `perturbation` sub-dict in the output JSON metadata.
- [x] **1.5** Removed `Get[PerturbativeReduction.wl]` from `_derive.py`.
- [x] **1.6** Deleted `_wls_phase_c_perturbative_reduction` function and its call site in `_derive.py`.
- [x] **1.7** Deleted `tidal/wolfram/PerturbativeReduction.wl` (525 lines, v5 JLM module).
- [x] **1.8** Deleted `tests/wolfram/test_perturbative_reduction.wl`.
- [x] **1.9** Rewrote `tests/test_phase_c_integration.py::TestWlsGeneration` — asserts metadata carries `"small_parameters" -> {...}` and `"perturbation_order" -> N` when `[perturbation]` configured; asserts no JLM artefacts remain.
- [ ] **1.10** Re-derive `torsion_gertsenshtein` and verify `order_in_eps` tags — deferred to Stage 6 (empirical validation block) to avoid derivation churn before Python side (Stages 2-5) can consume the tags.

### Smoke tests (2026-04-19)

- `wolframscript -file tests/wolfram/test_export_json.wls` → all 71 tests pass (exit 0).
- `ComputeOrderInEps` verified on synthetic inputs: `b5^2` → 2, `2*b5*eps + 3*b5^2` → 2, `5` → 0, empty params → 0.
- `BuildTermResult` verified: with `orderInEps=1`, result has `"order_in_eps" -> 1` key; with default (`0`), key is omitted.
- `uv run pytest tests/test_cli.py -x` → 175 passed.
- `uv run pytest tests/test_phase_c_integration.py -x` → 33 passed.
- Manual WLS inspection for `gertsenshtein` with `[perturbation]`: metadata carries `"small_parameters" -> {kappa}` and `"perturbation_order" -> 1`, no JLM residue.

---

## Stage 2: Python JSON-loader order_in_eps support

**Goal**: Propagate `order_in_eps` from JSON into `OperatorTerm` and provide filtering on `EquationSystem`.

**Estimated**: 1 day.

### Tasks

- [ ] **2.1** Add `order_in_eps: int = 0` field to `OperatorTerm` dataclass in `tidal/symbolic/json_loader.py`.
- [ ] **2.2** Parse `order_in_eps` from JSON in `OperatorTerm.from_dict` (default 0 for backward compat).
- [ ] **2.3** Add `EquationSystem.filter_by_order(n: int) -> EquationSystem` returning copy with only `order_in_eps == n` terms.
- [ ] **2.4** Add `EquationSystem.max_order() -> int` and `has_corrections() -> bool`.
- [ ] **2.5** Unit tests in `tests/test_json_loader.py` with synthetic mixed-order fixtures — verify round-trip, filter_by_order, has_corrections.

---

## Stage 3: Pass-0 eigendata export from modal solver

**Goal**: Modal solver returns eigenvalues/eigenvectors/Schur operators alongside `SolverResult` for use by Pass 1.

**Estimated**: 1 day.

### Tasks

- [ ] **3.1** Add `return_eigendata: bool = False` parameter to `solve_modal` in `tidal/solver/modal.py`.
- [ ] **3.2** When True, return `eigendata` dict with keys: `V`, `D_diag`, `V_inv`, `alpha = V_inv @ y0_hat`, `mode_k`, `schur_ops`, `state_layout`.
- [ ] **3.3** `schur_ops` must contain the S_cc⁻¹·S_cd matrices per mode computed during Pass 0 constraint elimination (modal.py:1020-1099).
- [ ] **3.4** Unit test `tests/test_modal_eigendata.py`: verify `V · diag(exp(D·t)) · V⁻¹ · alpha` reconstructs Pass 0 `SolverResult` at all timesteps to machine precision.

---

## Stage 4: Modal solver closed-form Duhamel — THE CRITICAL WORK

**Goal**: Pass 1 solves `dy/dt = A·y + M_src·y⁰` in closed form via φ₁ kernel. Constraint fields in source terms Schur-substituted into dynamical state.

**Estimated**: 5 days.

### Tasks

- [ ] **4.1** Add `_PHI1_DEGENERACY_THRESHOLD = 1e-5` module constant with Al-Mohy & Higham 2011 §3 Table 3.1 citation in docstring.
- [ ] **4.2** Implement `_duhamel_kernel(lam, mu, t)` with direct formula `exp(λt)·expm1(z)/z` for `|z| > threshold` and 12-term Taylor `t·exp(λt)·Σ z^k/(k+1)!` otherwise.
- [ ] **4.3** Implement `_build_source_matrix_k(correction_spec, k_vec, schur_ops, state_layout) -> NDArray`. For each term in each correction equation: if dynamical field → direct slot; if constraint field → expand via Schur row.
- [ ] **4.4** Implement `_evolve_per_mode_with_linear_source(M_src_k, y0_k, t_eval, eigendata_k)` — vectorised G-kernel evaluation per mode.
- [ ] **4.5** Extend `solve_modal` with `linear_source_matrices: dict[int, NDArray] | None = None` and `eigendata: dict | None = None` parameters. When provided: Duhamel path; when None: existing homogeneous path (unchanged).
- [ ] **4.6** Constraint-field correction recovery: after Pass 1 dynamical solve, recover constraint field correction via Schur_base + Schur_source augmentation.
- [ ] **4.7** Tests in `tests/test_modal_duhamel.py`:
  - Driven harmonic oscillator (analytical) to 1e-14
  - Constant source
  - Zero source fallback
  - Exact degeneracy `G(λ, λ; t) = t·exp(λt)`
  - Cross-field coupling
  - Constraint-field Schur substitution
- [ ] **4.8** Tests in `tests/test_modal_duhamel_degeneracy.py`: `|μ−λ|·t` sweep from 1e-15 to 1.0, relative error vs mpmath 50-digit reference ≤ 1e-13.

---

## Stage 5: Perturbative solver driver + CLI

**Goal**: `PerturbativeSolver` orchestrates Pass 0, Pass 1, combines. CLI flag `--perturbative-order N`.

**Estimated**: 5 days.

### Tasks

- [ ] **5.1** Create `tidal/solver/perturbative_driver.py` with `PerturbativeSolver` class and `PerturbativeResult` dataclass.
- [ ] **5.2** Implement `PerturbativeSolver.solve(ic, t_grid, order=1, params=None) -> PerturbativeResult`:
  - No corrections path: delegate to ModalSolver directly
  - Pass 0: `ModalSolver(base_spec).solve(..., return_eigendata=True)`
  - Pass n: build source matrices from `full_spec.filter_by_order(n)`, call `solve_with_source`
  - Combine via `total = Σ eps^n · q_n` evaluated at simulation params
- [ ] **5.3** Implement `_compute_validity(eigendata, full_spec, params)`:
  - `validity_param = max(ε_i · ω_max² · t_end)` across small params
  - Return dict with `omega_max`, `m_scalaron_min`, `validity_param`, `warn_level`
- [ ] **5.4** Add `--perturbative-order N` flag to `tidal simulate`. Default: 0 if `[perturbation]` absent; 1 if present. Reject `order > max_order_in_json`.
- [ ] **5.5** Wire `PerturbativeSolver` into `_simulate.py` dispatch. Validity warnings routed through `error_with_hint`.
- [ ] **5.6** Tests in `tests/test_perturbative_driver.py`: scalar toy (analytical match), synthetic 2-field with correction, validity warning firing/not firing.

---

## Stage 6: Validation + Ostrogradsky removal

**Goal**: End-to-end physics validation on toy theories and real examples. Delete `tidal/symbolic/ostrogradsky.py`.

**Estimated**: 4 days.

### Tasks

- [ ] **6.1** Algebraic-constraint toy theory test: base φ dynamical + h constraint (h = g·φ), correction b₅·(∂h)² + couplings. Verify h¹ ≠ 0, φ¹ matches analytical O(b₅) correction. `tests/test_perturbative_constraint_toy.py`.
- [ ] **6.2** Parker-Simon FLRW regression: 1+0D with a⁽⁰⁾(t) = t^(1/2), solve a⁽¹⁾ via Duhamel, match Parker-Simon Eq. 3.39-3.44 (gr-qc/9211002) to 1e-12. `tests/test_perturbative_parker_simon.py`.
- [ ] **6.3** Classical driven oscillator: `ÿ + 2γẏ + ω₀²y = A·exp(iω_d·t)`, sweep (γ, ω₀, ω_d) incl. near-degeneracy ω_d ≈ ω₀ − γ. Match to 1e-14. `tests/test_perturbative_driven_oscillator.py`.
- [ ] **6.4** Euler-Heisenberg refractive indices: verify `n²_⊥ − 1 = 4ρB²`, `n²_∥ − 1 = 7ρB²` emerge from Pass 1. `tests/test_perturbative_euler_heisenberg.py`.
- [ ] **6.5** R̃² PGT regression: h⁽¹⁾_{4,7,9} non-zero, energy conservation preserved (dE/E ≤ 1e-2), b₅=0 reproduces pre-perturbation baseline. `tests/test_perturbative_rtilde.py`.
- [ ] **6.6** Delete `tidal/symbolic/ostrogradsky.py` (520 lines).
- [ ] **6.7** Replace Ostrogradsky call in `json_loader.py:1239-1251` with recommendation warning when `time_order > 2 and not spec.has_corrections()`.
- [ ] **6.8** Remove `fields=` kwarg plumbing from `tidal/measurement/_energy.py` and `tidal/measurement/_conversion.py`.
- [ ] **6.9** Re-derive the three theories currently relying on Python Ostrogradsky: `graviton_torsion`, `torsion_gertsenshtein`, `torsion_gertsenshtein_combined`. Verify all load/simulate.
- [ ] **6.10** Close issues #195, #196. Move #164 anchor.

---

## Stage 7: Documentation

**Goal**: Update .tex docs + MEMORY.md. Version bump.

**Estimated**: 3 days.

### Tasks

- [ ] **7.1** Rewrite `docs/tex/perturbative_reduction.tex` to v6 iterative approach. Move JLM discussion to historical appendix.
- [ ] **7.2** Create `docs/tex/perturbative_numerical.tex` — iterative driver, modal Duhamel derivation, constraint-field Schur substitution, validity monitoring.
- [ ] **7.3** Update `docs/tex/architecture.tex` — Phase C is now order tagging + solver driver (not JLM).
- [ ] **7.4** Update MEMORY.md / `perturbative_reduction_research.md` — mark v6 complete.
- [ ] **7.5** CHANGELOG + minor version bump.
- [ ] **7.6** Final validation: full test suite + re-derive all 20 examples + baseline match.

---

## Verification Checklist

- [ ] `tidal/wolfram/PerturbativeReduction.wl` deleted
- [ ] `tests/wolfram/test_perturbative_reduction.wl` deleted
- [ ] ExportJSON.wl emits `order_in_eps` per OperatorTerm when `[perturbation]` configured
- [ ] `_derive.py` line 372 (Get[PerturbativeReduction.wl]) and line 4363 (`_wls_phase_c_perturbative_reduction`) deleted
- [ ] `EquationSystem.filter_by_order(n)`, `.max_order()`, `.has_corrections()` work on synthetic fixtures
- [ ] Pass-0 eigendata export includes `schur_ops` and reconstructs Pass 0 `SolverResult` to machine precision
- [ ] φ₁ degeneracy sweep: error ≤ 1e-13 from 1e-15 to 1.0 in `|μ−λ|·t`
- [ ] Constraint-field Schur substitution in M_src produces correct source matrix (toy theory test)
- [ ] Driven oscillator Duhamel matches analytical to 1e-14
- [ ] Parker-Simon FLRW regression matches to 1e-12
- [ ] Euler-Heisenberg `n²−1 = 4ρB²` emerges from Pass 1
- [ ] R̃² PGT: h⁽¹⁾_{4,7,9} non-zero, energy conservation preserved, b₅=0 baseline unchanged
- [ ] `tidal/symbolic/ostrogradsky.py` deleted; warn issued for time_order > 2 without perturbative flag
- [ ] `fields=` kwarg removed from `_energy.py` and `_conversion.py`
- [ ] All 20 example theories re-simulate at `--perturbative-order=0` matching baselines
- [ ] Full test suite passes
- [ ] Documentation updated in `docs/tex/`

---

## Timeline Summary

| Stage | Description | Days | Status |
|-------|-------------|------|--------|
| 0 | Post-E-L lagComp dual application | 2 | ✅ complete |
| 1 | Order-tagging in ExportJSON.wl + delete v5 JLM | 3 | ⏳ pending |
| 2 | Python JSON-loader order_in_eps support | 1 | ⏳ pending |
| 3 | Pass-0 eigendata export (incl. schur_ops) | 1 | ⏳ pending |
| 4 | Modal solver closed-form Duhamel | 5 | ⏳ pending |
| 5 | Perturbative solver driver + CLI flag | 5 | ⏳ pending |
| 6 | Validation + Ostrogradsky removal | 4 | ⏳ pending |
| 7 | Documentation | 3 | ⏳ pending |

**Remaining: ~22 days (~4.5 weeks)**

---

## Historical Note — v5 (JLM) Deprecated

Prior Stages 1-3 under the JLM approach (symbolic order-reduction via substitution of base EOM into correction terms) were implemented but are being replaced. Reasons:
- JLM cannot handle algebraic-constraint fields promoted to dynamical by correction (R̃² PGT h_4/h_7/h_9)
- Required theory-specific classification (gauge fixing, scalarization rewriter)
- Secular growth error bound not cleanly stated

v6 iterative-numerical approach:
- Ghost-free by construction (LHS operator always 2nd-order base)
- Theory-agnostic (no classification)
- Machine-precision closed-form Duhamel via φ₁ kernel
- Sharp error bound `C·ε²·ω²·t` for linear theories (beats FKY's nonlinear `C·ε²·(t·M_UV)³`)

Supersession date: 2026-04-19. PerturbativeReduction.wl to be deleted in Stage 1.
