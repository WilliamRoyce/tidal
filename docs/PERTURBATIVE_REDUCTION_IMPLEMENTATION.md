# Perturbative Reduction Implementation Tracking (v6)

**Status (2026-04-19): COMPLETE.** All 19 sub-tasks across 6 phases done (6C.2 Euler-Heisenberg example deferred to issue #271 — xAct `Validate::repeated` on quartic `(F·F)²` contraction). Full Python suite: **1920 passed, 66 skipped**. v0.33.0 release.

Commit trail (branch `feature/torsion-sweeps`):

| Phase | Commit | Summary |
|-------|--------|---------|
| 6A.1 | `d30c42f` | Gap B — AST-based LHS demotion for base_spec |
| 6A.2 | `65cac69` | Gap C — Pass 1 constraint-field Schur recovery |
| 6A.3 | `54aab3a` | Algebraic-constraint toy theory (validates B+C) |
| 6B.1+6B.2 | `ceb6e63` | Delete ostrogradsky.py + migration error + fields= kwarg |
| 6B.3 | `689b182` | Re-derive torsion_gertsenshtein with [perturbation] |
| 6B.4 (part 1) | `86f9087` | Re-derive graviton_torsion |
| 6B.4 (part 2) | `7ab47a5` | Re-derive torsion_gertsenshtein_combined |
| 6C.1 | `e567157` | R̃² PGT regression test |
| 6C.3+6C.4 | `2e3d297` | Near-degeneracy + cross-block + review cleanup |
| Stage 7 | _(in progress)_ | Docs, CHANGELOG, version bump |

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

## Stage 2: Python JSON-loader order_in_eps support — COMPLETE

**Goal**: Propagate `order_in_eps` from JSON into `OperatorTerm` and provide filtering on `EquationSystem`.

**Estimated**: 1 day. **Actual**: <1 day (2026-04-19).

### Tasks

- [x] **2.1** Added `order_in_eps: int = 0` field to `OperatorTerm` dataclass in `tidal/symbolic/json_loader.py:219`. Documented in the docstring.
- [x] **2.2** `OperatorTerm.from_dict` parses `order_in_eps` with `.get("order_in_eps", 0)` default for backward compat with pre-v6 JSONs.
- [x] **2.3** `EquationSystem.filter_by_order(n)` at `json_loader.py:1144`: returns a copy of the system with each equation's `rhs_terms` filtered to those matching `order_in_eps == n`. Preserves equations and LHS structure so the state layout is unaffected.
- [x] **2.4** Added `EquationSystem.max_order() -> int` and `has_corrections() -> bool`.
- [x] **2.5** Tests: `TestOperatorTerm::test_order_in_eps_{default_zero,parsed_from_json,second_order,direct_construction}` and `TestEquationSystemOrderInEps::{test_has_corrections_baseline_false, test_has_corrections_true_with_mixed_orders, test_filter_by_order_zero_keeps_base_terms_only, test_filter_by_order_one_keeps_correction_terms_only, test_filter_by_order_preserves_immutability, test_filter_by_order_empty_high_order}`. Synthetic `_mixed_order_spec_data()` fixture mirrors the R̃² PGT structure (base + b₅ correction source).
- Side-fix: updated `tests/test_wls_helpers.py::TestWlZeroComponent` to match Stage 0 dual-application (lagComp block added to output).

### Smoke tests (2026-04-19)

- `uv run pytest tests/test_json_loader.py -x -q` → 105 passed.
- `uv run pytest tests/ -x -q --ignore=tests/integration` → 1811 passed, 66 skipped (full Python suite, no regressions).

---

## Stage 3: Pass-0 eigendata export from modal solver — COMPLETE

**Goal**: Modal solver returns eigenvalues/eigenvectors/Schur operators alongside `SolverResult` for use by Pass 1.

**Estimated**: 1 day. **Actual**: <1 day (2026-04-19).

### Tasks

- [x] **3.1** Added `return_eigendata: bool = False` keyword to `solve_modal`. Threaded through a new `collect_eigendata: list[dict] | None = None` parameter on `_evolve_per_mode` — a mutable container the evolution routine appends each block's eigendecomposition into, avoiding a breaking change to the evolution function's return signature.
- [x] **3.2** When `return_eigendata=True`, the returned `SolverResult` carries an `"eigendata"` dict with keys: `blocks` (list of per-block dicts with `slot_indices`, `V`, `D_diag`, `V_inv`, `alpha = V_inv @ y0_hat`), `mode_k` (the k-grid), `state_layout` (the *reduced* layout when Schur eliminate runs; full layout otherwise).
- [x] **3.3** `schur_ops` sub-dict emitted when the constraint-field path runs: `recovery_matrix` (the `-S_cc⁻¹·S_cd` per mode, already computed at modal.py:1020-1099), `constraint_field_names` tuple, and `orig_to_reduced` slot-index map. Pass 1 (Stage 4) uses this to Schur-substitute constraint-field references when assembling `M_src(k)`.
- [x] **3.4** Six-test suite at `tests/test_solver_modal.py::TestEigendataExport`:
  - `test_eigendata_key_absent_by_default` — baseline no-eigendata behaviour
  - `test_eigendata_structure_present_when_requested` — shape of the dict
  - `test_eigendata_invertibility` — `V @ V_inv ≈ I` per mode
  - `test_eigendata_alpha_matches_vinv_y0` — α derivation correct
  - `test_eigendata_reconstructs_pass_zero` — `y(t) = V·diag(exp(D·t))·α` matches snapshot to < 1e-12
  - `test_position_dependent_raises` — `NotImplementedError` on Krylov path

### Smoke tests (2026-04-19)

- `uv run pytest tests/test_solver_modal.py -q` → 44 passed (38 pre-existing + 6 new).
- `uv run pytest tests/ -x -q --ignore=tests/integration` → 1817 passed, 66 skipped (full Python suite, no regressions from the refactor).
- Position-dependent path correctly rejects `return_eigendata=True` with a clear error pointing at Pass 1 Duhamel's constant-coefficient requirement.

---

## Stage 4: Modal solver closed-form Duhamel — COMPLETE

**Goal**: Pass 1 solves `dy/dt = A·y + M_src·y⁰` in closed form via φ₁ kernel. Constraint fields in source terms Schur-substituted into dynamical state.

**Estimated**: 5 days. **Actual**: 1 day (2026-04-19).

### Tasks

- [x] **4.1** Added `_PHI1_DEGENERACY_THRESHOLD = 1e-5` in `tidal/solver/modal.py` with Al-Mohy & Higham 2011 §3 Table 3.1 citation in the module-level docstring block that precedes the kernel.
- [x] **4.2** Implemented `_duhamel_kernel(lam, mu, t)` with the correct direct formula `exp(λt) · t · expm1(z) / z` (fixed mid-review — initial version was off by a factor of t) and 12-term Taylor fallback `t · exp(λt) · Σ z^k / (k+1)!` for `|z| ≤ threshold`. Handles scalar and broadcast-compatible array inputs.
- [x] **4.3** Implemented `_build_source_matrix_k(correction_spec, layout, coeff_eval, k_grid, rfft_shape, schur_ops=None)`. For each RHS term in the correction spec: dynamical target → direct slot write; constraint target → expanded via `recovery_matrix[m, c_idx, :]` row from the Pass 0 Schur eliminate. Row selection mirrors `_build_per_mode_matrices` (velocity slot for 2nd-order equations, field slot for 1st-order).
- [x] **4.4** Implemented `_evolve_duhamel_per_mode(eigendata, M_src_k, t_eval, layout, grid)`. Block-aware: for each Pass 0 block, projects `M_src` sub-block into eigenbasis (`β = V⁻¹·M_src·V`), evaluates `z_i(t) = Σ_j β_ij · α_j · G(λ_i, λ_j; t)` via the Duhamel kernel, transforms back to slot space, and inverse-FFTs to physical space. Guards against cross-block coupling (raises `NotImplementedError` if the correction spec links previously-independent sectors).
- [x] **4.5** Added top-level `solve_modal_pass1(eigendata, correction_spec, grid, t_eval, *, parameters=None) -> SolverResult`. Thin wrapper over `_build_source_matrix_k` + `_evolve_duhamel_per_mode`. Independent of `solve_modal` so the driver (Stage 5) can compose the two cleanly.
- [ ] **4.6** Constraint-field correction recovery: deferred to Stage 5. The Pass 1 state currently returns only the dynamical (reduced) layout; the Stage 5 driver will apply `recovery_matrix @ y_dyn¹` (base Schur recovery) plus the `Schur_source_correction` augmentation for constraint fields like h_4 in R̃² PGT. No correction-recovery logic exists yet because the KG toy has no constraints.
- [x] **4.7** `tests/test_modal_duhamel.py::TestSolveModalPass1AnalyticalMatch` (3 tests):
  - Driven Klein-Gordon with ε·φ correction: q⁽⁰⁾ + q⁽¹⁾ matches the full-ε solver result to O(ε²) ≈ 1e-2 at ε=0.05, t=2. Also asserts a >4× improvement over bare Pass 0.
  - Empty correction → zero Pass 1 output.
  - IC is zero: `y⁽¹⁾(t_eval[0]) = 0` per construction.
- [x] **4.8** `tests/test_modal_duhamel_degeneracy.py` (40 tests): scalar edge cases (exact degeneracy, t=0, broadcast vectorized) + parametrised sweep with `log10|μ−λ|` from -15 to 0 on real and imaginary axes, t ∈ {0.1, 1.0, 10.0}, plus crossover-continuity check. mpmath 50-digit reference; all pass with relative error ≤ 1e-12.

### Smoke tests (2026-04-19)

- `uv run pytest tests/test_modal_duhamel_degeneracy.py -q` → 40 passed.
- `uv run pytest tests/test_modal_duhamel.py -q` → 3 passed.
- `uv run pytest tests/ -x -q --ignore=tests/integration` → **1860 passed**, 66 skipped (no regressions from the prior 1817 baseline).

### Deferred for Stage 5

1. Constraint-field correction recovery (Task 4.6). Requires the driver pipeline to compose Pass 0 reconstruction, Pass 1 dynamical evolution, and the augmented Schur recovery for constraint fields.
2. Full R̃² PGT regression test. The toy KG test validates the kernel path; the real PGT run needs the driver (Stage 5) + re-derived JSON with order tags.

### New dependency

- `mpmath` added to dev dependencies (`pyproject.toml` via `uv add --dev mpmath`) — used by the degeneracy sweep as a 50-digit reference for the Duhamel kernel. Not imported by production code.

---

## Stage 5: Perturbative solver driver + CLI — COMPLETE

**Goal**: `PerturbativeSolver` orchestrates Pass 0, Pass 1, combines. CLI flag `--perturbative-order N`.

**Estimated**: 5 days. **Actual**: 1 day (2026-04-19).

### Tasks

- [x] **5.1** Created `tidal/solver/perturbative_driver.py` with `PerturbativeResult` dataclass (per-order SolverResults, combined `total`, `validity` dict) and `PerturbativeSolver` class.
- [x] **5.2** `PerturbativeSolver.solve(y0, grid, t_span, *, order, parameters, num_snapshots, small_parameters)` implemented. Pass 0 via `solve_modal(..., return_eigendata=True)`; Pass n via `solve_modal_pass1` for each order. No-correction path delegates cleanly; cross-block coupling detected in the Duhamel routine.
- [x] **5.3** `_compute_validity(eigendata, small_parameters, parameters, t_end)` returns `{omega_max, eps_values, validity_param, dominant_parameter, warn_level}` with thresholds `_VALIDITY_WARN_THRESHOLD = 0.1` and `_VALIDITY_ERROR_THRESHOLD = 1.0`, referenced to the plan's "TIDAL's parameter space" table.
- [x] **5.4** `--perturbative-order N` added to `tidal simulate`. Default rule: `1` when JSON metadata has a `perturbation` block; `0` otherwise. Explicit values override the default.
- [x] **5.5** `tidal/cli/_simulate.py` modal branch splits: when `pert_order > 0 and spec.has_corrections()` it routes through `PerturbativeSolver`; otherwise the existing `solve_modal` path runs unchanged. Validity warnings surface via `_cwarn` in two bands (warn @ 0.1, error @ 1.0).
- [x] **5.6** Tests:
  - `tests/test_perturbative_driver.py::TestPerturbativeSolverAPI` — `max_order`, `has_corrections`, rejects over-order.
  - `TestPerturbativeSolverSolve` — order=0 returns base-only; order=1 adds non-zero correction; validity monitor flags large ε·ω²·t at the error band.
  - `TestPerturbativeSolverValidity` — small-ε stays in "ok"; parameter map populated correctly.
  - `TestPerturbativeCLIFlag` — CLI invocations exercise explicit order=1, default order=1 (from metadata), and order=0 (bypass). 11 tests total.
- Also exported `PerturbativeSolver`, `PerturbativeResult`, `solve_modal`, `solve_modal_pass1` from `tidal.solver`.

### Smoke tests (2026-04-19)

- `uv run pytest tests/test_perturbative_driver.py -q` → 11 passed.
- `uv run pytest tests/ -q --ignore=tests/integration` → **1871 passed**, 66 skipped (no regressions from the 1860 baseline).
- `uv run ruff check ...` on all Stage 5 files → clean.

### Deferred to Stage 6

- Full constraint-field correction recovery (augmented Schur step including source contribution for `h_4¹` etc.). Current driver zero-fills constraint slots in the Pass 1 output, which is correct when the correction source does not back-propagate into the constraint equations (true for the KG toy, not for R̃² PGT). Stage 6 adds the augmented recovery.
- Real R̃² PGT re-derivation + regression test against the base b₅=0 baseline.

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

- [x] **7.1** Rewrite `docs/tex/perturbative_reduction.tex` to v6 iterative approach. Move JLM discussion to historical appendix. (Completed 2026-04-19; extended 2026-04-20 with §Power-of-Contraction Normalisation and Matter-Only Derivative Dependence for issue #271.)
- [x] **7.2** Create `docs/tex/perturbative_reduction_design.tex` — iterative driver, modal Duhamel derivation, constraint-field Schur substitution, validity monitoring. (Filename finalised as `perturbative_reduction_design.tex`; engineer-facing companion to §7.1. Both files wired into `main.tex` 2026-04-20.)
- [x] **7.3** Update `docs/tex/architecture.tex` — Phase C is now order tagging + solver driver (not JLM).
- [x] **7.4** Update MEMORY.md / `perturbative_reduction_research.md` — mark v6 complete.
- [x] **7.5** CHANGELOG + minor version bump. (Version line reached v0.33.9 at session close.)
- [x] **7.6** Final validation: full test suite + re-derive representative examples + baseline match. (1 991 Python tests green; `gertsenshtein` $\approx 49$~s, `coupled_scalars` $\approx 14$~s, `euler_heisenberg` $\approx 19$~s (newly derivable); only `derivation_hash` fields change on existing JSONs.)

### Stage 7 Closure Entry — Issue #271 (2026-04-20, v0.33.9)

**Problem.** Euler–Heisenberg $(F_{\mu\nu}F^{\mu\nu})^2$ refused to derive. Investigation during session identified two orthogonal root causes.

**Root cause 1 — Power auto-distribution.** Mathematica's built-in rule `Power[Times[a,b,...], n] → a^n·b^n·...` fires at assignment time on the scalar contraction, producing repeated abstract indices ⇒ xAct `Validate::repeated + Throw[Null]`. Fix in `tidal/cli/_derive.py::_wls_lagrangian` (lines 986–1017): wrap user Lagrangian in `Hold[]`, apply ReplaceRepeated `Power[X, n] → Scalar[X]^n` (where `X` carries abstract indices and is not already a `Scalar`), then `ReleaseHold`. xAct treats `Scalar[...]^n` as opaque — distribution blocked, `DecomposeScalarExpression` visits each instance with a fresh index scope. Commit 9d9e73f.

**Root cause 2 — Matter-only CD precompute gate.** `_wls_precompute_cd_component_values` skipped the CD ComponentValue precomputation for `len(dyn_fields) < 2`, a performance heuristic. For EH (single dynamical photon appearing only via `F = dA`) this left `CD1eqca[{i,-cart},{j,-cart}]` unresolved in `lagComp`, and Component E-L found zero field functions. Fix in `tidal/cli/_derive.py` (lines 1513–1620): correctness-aware gate runs precompute whenever `len(dyn_fields) ≥ 2` OR any derived_field with `"CD["` in its definition is used OR any dyn-field/matter-pert source field appears as `CD[...][name[...]]` in the user expression. Commit 830a442.

**Supplementary fix — `$CDShorthandReverseRules` safety net.** These rules (generated at `_derive.py:1407` since v0.28) were never applied. Added application to `lagComp` before Component-E-L field-function detection (lines 5857–5869), FreeQ-gated to keep it O(LeafCount) when the precompute already resolved everything.

**Documentation.** §Power-of-Contraction Normalisation and Matter-Only Derivative Dependence in [`docs/tex/perturbative_reduction.tex`](../docs/tex/perturbative_reduction.tex); engineer-facing counterpart §Power-of-Contraction Normalisation — Implementation Details in [`docs/tex/perturbative_reduction_design.tex`](../docs/tex/perturbative_reduction_design.tex); timing row in [`docs/tex/derivation_performance.tex`](../docs/tex/derivation_performance.tex) §Post-v6 Reference Timings.

**Verification.**

| Theory | Wall-time | Components | Notes |
|---|---|---|---|
| `euler_heisenberg` | ~19s | 4 | newly supported; `order_in_eps=1` with `2·B0²·ρ`, `-2·B0²·ρ`, `-6·B0²·ρ` |
| `gertsenshtein` | ~49s | 6 | baseline 46s, +7% within variance (FreeQ-gated safety net) |
| `coupled_scalars` | ~14s | 2 | unchanged |
| `scalar_potential_well` | ~10s | 1 | fast-path precompute skip preserved |

1,991 Python tests green. Existing JSON outputs change only in `derivation_hash`.

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

---

## Post-Audit Remediation (v0.33.x, 2026-04-19)

After v0.33.0 shipped, a two-pass critical audit identified 14 issues
(#272–#286, plus the pre-existing xAct blocker #271) ranging from silent-drop
correctness bugs (#272, #273) to validity-monitor science errors (#284, #285).
The remediation plan (`.claude/plans/flickering-gathering-orbit.md`, approved
2026-04-19) landed across R1–R7 commits:

- **R1 (correctness)** e4c1d1c: `_build_source_matrix_k` drops logged and
  surfaced on `PerturbativeResult.validity["correction_drops"]` (#272);
  `order >= 2` raises `NotImplementedError` with #273 reference;
  `base_spec`/`filter_by_order`/`normalize_kinetic_coefficients` recompute
  mass/coupling matrices (#274); cross-block threshold is scale-relative
  `max(1e-14, max|M|·1e-10)` (#275).
- **R2 (typed contracts)** bd8f3ae: `PerturbativePass1Result` TypedDict
  replaces `# type: ignore` smuggling (#279); `PerturbativeResult` carries
  `spec` + `full_spec`, CLI reads `pert_result.spec` instead of rebinding
  (#276).
- **R3 (runtime invariants)** bd8f3ae: silent `except ValueError` in
  `_resolve_scheme` dropped (#277); `solve_modal` rejects `time_order > 2`
  with actionable error (#278).
- **R4 (validity monitor)** bd8f3ae: formula corrected to `ε·ω·t` (#284);
  base-theory tachyon check added as separate `base_level` diagnostic
  (#285).
- **R5 (physics tests)** a6e5dd5: R̃² regression rewritten around non-gauge
  channels + drop-category assertions (#280); new baseline file with
  structural + synthetic-spec checks (#281, trajectory-level regression
  deferred to #289 pending derivation of a proper b5=0 reference JSON);
  Taylor-branch crossover pipeline-level test against mpmath 50-dps
  (#282); constraint-IC projection + JSON round-trip variants of the
  trivially-true / unreachable-path tests (#283).
- **R7 (rendering)** /tmp/fix_issues.py: unicode batch-fix of issue
  titles and bodies (ε / ω / φ / λ etc. now render properly in #272–#286).

### Phase 6B.5 retrospective (R6.1 / #286)

Retrospective check that the Ostrogradsky-removal commit (`ceb6e63`) did not
silently loosen live test coverage:

- Only two test files were modified in `ceb6e63`:
  - `tests/test_json_loader.py`: **added** 42 lines — new
    `test_high_order_without_perturbation_raises` verifying the migration
    error fires on synthetic `time_order=4` JSON without `[perturbation]`.
    No existing assertions touched.
  - `tests/test_latex.py`: **added** 25 lines — three theories marked
    `@pytest.mark.xfail(strict=True)` pending Phase 6B.3/6B.4 re-derivation.
    `strict=True` means XPASS becomes FAILED, so the marker would have
    auto-flipped the moment the migrations landed. Post-migration, the
    `_PENDING_V6_REDERIVATION` set has been emptied (the three theories
    now load cleanly through the post-Ostrogradsky loader). Verified by
    running `test_latex.py::TestAllExamplesRender::test_no_crash` —
    25/25 pass, 0 xfailed.
- **No assertions loosened**; no tests deleted outside of the
  `tidal/symbolic/ostrogradsky.py` deletion itself (which carried no
  standalone test file).
- **Hamiltonian consistency check (#178 still wired)**: all three
  re-derived JSONs (`torsion_gertsenshtein`, `graviton_torsion`,
  `torsion_gertsenshtein_combined`) carry `canonical.hamiltonian_terms`
  (78, 256, 237 entries respectively). The constrained-Hamiltonian logic
  referenced in `json_loader.py:1106` (`#178`) still fires at load time;
  the shipped JSONs pass the consistency check without the legacy
  `fields=` kwarg that was removed alongside Ostrogradsky in `ceb6e63`.
- **Pre-v6 vs post-v6 measurement regression**: not executed in-session
  because no pre-v6 archived `tidal measure` output exists in the repo
  to compare against. Documented as a backlog item rather than blocking
  — the theory-level invariants (all three JSONs load, run, and respect
  the iterative-perturbative split) plus the structural R5.2 baseline
  checks cover the migration correctness without needing a byte-for-byte
  pre/post measurement comparison.

### Follow-up issues carried forward

- **#271** — Euler-Heisenberg xAct blocker (`Validate::repeated` on
  `(F·F)²`). Out of scope for R1–R7 remediation; Wolfram-pipeline
  investigation.
- **#273** — order ≥ 2 Pass 2 Duhamel against `q¹` is now gated with
  `NotImplementedError`; proper implementation needs a new IC-amplitude
  contract.
- **#289** — Derive `torsion_gertsenshtein_b5_zero.json` so R5.2's
  trajectory-level baseline regression can be activated (the shipped
  `torsion_gertsenshtein_minimal_propagating.json` uses a different
  Lagrangian; not a b5=0 variant).
