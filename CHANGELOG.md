# Changelog

All notable changes to TIDAL are documented in this file. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
uses semantic versioning as interpreted by `scripts/bump_version.py`.

This changelog was introduced in v0.30.2 and backfilled for the
canonicalization-architecture trail (v0.29.3 → v0.30.1) because those
commits reworked a load-bearing part of the Wolfram pipeline that future
maintainers need to be able to trace. Earlier versions are not
retroactively covered; see `git log` for the full history.

## [v0.34.0 – v0.47.9] — not individually recorded

Changelog maintenance lapsed during the PGT survey campaign and the MSci
write-up (April–August 2026). The work in this range covers the v3 inference
architecture, the coupling-space survey and its HPC tooling, the localised-field
and cubed-sphere campaigns, solver stability guards, and the manuscript, poster
and talk build infrastructure. `git log v0.33.4..v0.47.9` is the authoritative
record for this window; entries below resume the per-version trail.

## [v0.33.4] — Pass 2 math documented (#273)

### Changed

- **Upgraded `order >= 2` gate diagnostic (#273)**. The
  `NotImplementedError` message now cites both the direct
  (`M_src_2 · y⁰(t)`, reuses Pass 1 machinery via
  `filter_by_order(2)`) and indirect (`M_src_1 · y¹(t)`, needs the
  triple-eigenvalue nested kernel `K_2(λ_i, λ_j, λ_k; t) =
  ∫₀ᵗ exp(λ_i(t-τ))·G(λ_j, λ_k; τ)dτ`) contributions.
  Analytical derivation plus degenerate-case handling is documented
  in `docs/tex/perturbative_reduction.tex` §Pass 2 and higher
  orders. The gate is preferred over a partial implementation
  because no shipped theory emits `order_in_eps=2` terms; the
  untested code path would silently produce wrong results.

### Added

- `docs/tex/perturbative_reduction.tex` §"Operator time-derivative
  scaling": documents the `λⁿ` factor applied per-time-order in
  `_evolve_duhamel_per_mode` (fixing #293). §"Pass 2 and higher
  orders" derives the analytical triple-eigenvalue kernel for
  future implementation.

## [0.33.3] — 2026-04-20 time-order λⁿ scaling + b5=0 reference JSON

### Fixed

- **`_build_source_matrix_k` dropped operator time_order (#293)**.
  Pass 1 correction terms with `d^n_t` operators targeting a
  dynamical field silently lost the `λⁿ` eigenvalue factor because
  `_EXACT_MULTIPLIERS` kept only the spatial multiplier. Latent
  bug — no shipped theory has a b5-coupled d^n_t term on a
  dynamical row (all route through the constraint augmented
  recovery path where `λⁿ` was already handled). Fix:
  - `_build_source_matrix_k` now returns `dict[int, ndarray]` keyed
    by operator time_order.
  - `_evolve_duhamel_per_mode` scales α by `λⁿ` per order before
    the Duhamel kernel.
  - New regression `TestPass1TimeDerivativeTargetingDynamical` in
    `tests/test_modal_duhamel.py`: damped KG at O(γ¹) matches
    `Φ = cos(ω₀t) + γ·[−t/2·cos(ω₀t) + 1/(2ω₀)·sin(ω₀t)] + O(γ²)`.
    Pre-fix this would fail by O(1).

### Added

- **b5=0 reference JSON (#289)** for the R5.2 trajectory-level
  regression. `examples/torsion_gertsenshtein/theory_b5_zero.toml`
  is a byte-level b5 → 0 reduction of `theory.toml` (removes the
  `b5·R̃²` Lagrangian term, the `b5` constant, and the
  `[perturbation]` section). Derived to
  `examples/data/torsion_gertsenshtein_b5_zero.json`: 38 fields, 0
  `order_in_eps` tags, `h_4`/`h_7`/`h_9` natively algebraic.
  New test `test_pass0_matches_b5_zero_reference` confirms the
  full theory's Pass 0 at b5=0 reproduces this reference
  field-for-field to rel_err < 1e-6.

## [0.33.2] — 2026-04-20 R8 augmented Schur constraint recovery

### Fixed

- **Pass 1 constraint recovery missing O(ε¹) contributions (#290)**.
  The flagship R̃² PGT theory (`torsion_gertsenshtein.json`) produced
  identically-zero Pass 1 output on `h_{4,7,9}` pre-R8 because all
  159 O(ε¹) corrections live on constraint-field rows and were
  silently dropped. The v6 plan explicitly anticipated this via
  Stage 4.6's `Schur_source_correction` augmentation ("deferred to
  Stage 5"), but the code was never written. R8 implements the
  augmented recovery:
  ```
  h_c¹ = recovery · y_dyn¹ + S_cc⁻¹ · [K · d^n_t(h_c⁰) − corr_1(y⁰)]
  ```
  - `tidal/solver/modal.py`: `_build_evolution_matrices` returns
    `S_cc_inv` and `singular_mask` alongside `recovery_matrix`;
    `solve_modal` attaches both to `eigendata["schur_ops"]`.
  - `tidal/symbolic/_kinetic_eval.py`: new `evaluate_at_one` and
    `evaluate_with_substitutions`; driver resolves kinetics at
    runtime so the ε factor is correctly baked in.
  - `tidal/solver/perturbative_driver.py`: new `_pre_demote_info`
    and `_compute_constraint_source_hat` compute the augmented
    source from Pass 0 eigendata. `_assemble_full_state_pass_n`
    applies `S_cc_inv @ source_hat` to the constraint slots.
  - Drop-record label for constraint-row corrections updated from
    `"row-missing"` to `"row-routed-to-augmented"` (no longer
    actually dropped).
- **Mpmath cross-check validation (R8.6)**. New
  `TestAugmentedRecoveryMatchesMpmath` with 3 tests validating
  the driver's h/phi ratio against the 2×2 dispersion light-mode
  eigenvector at 50 dps. At ε=1e-3, O(ε¹) truncation matches to
  rel_err < 1e-6; exact eigenvector to < 1e-5 (the O(ε²) gap).
  At ε=1e-2, truncation matches to < 1e-4.
- **Synthetic spec's `|h − g·phi|` now matches analytical
  `−ε·g·ω²·phi⁰`** to rel_err 1.8e-14 (machine precision). Pre-R8
  this was 2e-16, 14 orders of magnitude below the physics.

### Tracked

- **#273** order ≥ 2 Pass 2 Duhamel — unchanged; augmented
  recovery at O(ε¹) is a prerequisite but Pass 2 against q¹(t)
  remains a separate architectural extension.
- **#289** b5=0 reference JSON — unchanged; activates the R5.2
  trajectory-level baseline test.
- **#293** latent time_order bug in `_build_source_matrix_k` —
  uncovered during R8.7 investigation; doesn't affect shipped
  theories (all O(ε¹) corrections route through the constraint-
  row augmented path, where time_order IS handled correctly).

---

## [0.33.1] — 2026-04-19 post-audit remediation

### Fixed

- **Silent drop of correction terms (#272)**. `_build_source_matrix_k`
  now returns a drops list; unroutable correction terms surface as
  structured records on `PerturbativeResult.validity["correction_drops"]`
  and log a WARN per drop. The R̃² PGT regression test asserts every
  drop is of category `"row-missing"` (legitimate O(ε²) augmentation);
  any `"column-missing"` drop is a bug.
- **order ≥ 2 silently reused Pass 0 eigendata (#273)**. Now raises
  `NotImplementedError` with a reference to the tracking issue. The
  docstring no longer claims Pass 2 is API-accepted.
- **Stale mass/coupling matrices after `base_spec` (#274)**. `base_spec`,
  `filter_by_order`, and `normalize_kinetic_coefficients` now recompute
  the matrices from their filtered equations. `EquationSystem.__post_init__`
  no longer fires a UserWarning on every call.
- **Cross-block threshold false-triggers on ill-conditioned Schur (#275)**.
  The guard uses a scale-relative threshold
  `max(1e-14, max|M_src| · 1e-10)` instead of absolute 1e-14.
- **CLI spec-rebind workaround (#276)**. `PerturbativeResult` carries
  `spec` (base_spec) and `full_spec` attributes; CLI reads
  `pert_result.spec` instead of rebinding.
- **Silent `except ValueError` in `_resolve_scheme` (#277)**. Removed.
  `base_spec`'s actionable error now propagates to the user instead
  of being masked by modal-ineligibility fallback.
- **`solve_modal` time_order > 2 silent acceptance (#278)**. Rejected
  with actionable error pointing at `PerturbativeSolver` + the docs.
- **`PerturbativePass1Result` TypedDict (#279)**. Replaces the `# type:
  ignore` smuggling of `y_hat_dyn` and `correction_drops` into
  `SolverResult`. Driver's defensive `None`-skip on `y_hat_dyn` becomes
  an explicit `assert` (required field, not optional).
- **Validity monitor formula `ε·ω²·t` → `ε·ω·t` (#284)**. The secular
  phase-error bound. The old formula was conservative at ω > 1 but
  dangerously permissive at ω < 1.
- **Validity monitor blind to base-theory tachyons (#285)**. Adds a
  separate `base_level` / `base_stability_param` diagnostic (from
  `max Re(λ) · t_end`); the overall `warn_level` is now the worse of
  correction-side and base-stability-side bands.

### Changed

- **R̃² PGT regression rewritten around physics (#280)**. The prior
  test only asserted `isfinite + len==2 + warn=='ok'`. Replaced with
  non-gauge-channel assertions: photon IC excites Gertsenshtein
  h-constraint mixing, drops are all row-missing (no column-missing
  bugs), Pass 1 is identically zero on dynamical fields (theory-level
  O(ε) invariant of the shipped JSON — all b5 couplings live on
  constraint rows), Pass 1 amplitude scales linearly in b5 where
  non-zero, Pass 0 is b5-independent per-field.
- **Taylor-branch crossover pipeline-level test (#282)**. The previous
  `TestPass1NearDegeneracy` tested resonance, not near-degeneracy.
  Renamed to `TestPass1Resonance`; new `TestPass1NearDegeneracy`
  sweeps `|Δλ|·t` across {2e-11, 1e-5, 1e-2, 1} and asserts agreement
  with an mpmath 50-dps reference to rel_err < 1e-10 through
  `_evolve_duhamel_per_mode`.

### Added

- **Structural `base_spec` baseline checks (#281)**. New
  `tests/test_perturbative_rtilde_baseline.py` validates the full
  theory's `base_spec(['b5'])` against structural invariants
  (no order > 0 RHS terms, no b5 in any `coefficient_symbolic`,
  time_orders ≤ 2, idempotent under re-application, Pass 0 trajectory
  b5-independent across every common dynamical field). Trajectory-
  level regression against a true b5=0 reference JSON tracked in #289
  (the shipped `minimal_propagating` JSON is a different Lagrangian,
  not a b5→0 limit).
- **Constraint-IC projection test + JSON round-trip post-check test
  (#283)**. Replaces the prior trivially-true and unreachable-path
  unit tests with realistic user-facing flows.

### Tracked

- **#271** Euler-Heisenberg xAct `Validate::repeated` blocker (Wolfram
  pipeline; out of scope for R1–R7).
- **#289** Derive `torsion_gertsenshtein_b5_zero.json` so R5.2's
  trajectory-level baseline regression can be activated.

## [0.33.0] — 2026-04-19

### Added

- **v6 iterative perturbative reduction** (issue #267). Higher-
  derivative theories declaring `[perturbation] small_parameters=[...]`
  in their `theory.toml` are now evolved by a Parker-Simon iterative
  scheme instead of the legacy mechanical Ostrogradsky reduction.
  - `PerturbativeSolver` class in `tidal/solver/perturbative_driver.py`
    orchestrates Pass 0 (base, ε=0) + Pass 1 (closed-form Duhamel
    correction) using shared eigendata.
  - `tidal simulate --perturbative-order N` flag. Defaults to 1 when
    the JSON metadata carries a `perturbation` block, 0 otherwise.
  - `EquationSystem.base_spec(small_parameters)` demotes LHS kinetic
    coefficients that vanish at ε=0 to algebraic constraints, giving
    a ghost-free 2nd-order base system (Gap B).
  - `tidal/solver/modal.py::solve_modal_pass1` evaluates the closed-
    form φ₁-kernel Duhamel integral with Al-Mohy-Higham (2011) Taylor
    fallback near degenerate eigenvalues; validated against 50-digit
    mpmath reference to relative error < 1e-13 across |μ−λ|·t ∈
    [1e-15, 1].
  - Constraint fields recovered at Pass 1 via the existing Schur
    operator applied to the Pass 1 dynamical Fourier output
    (`c_hat = recovery_matrix @ y_hat_dyn`) — Gap C.
  - Validity monitor: flags `ε · ω² · t_end > 0.1` as "warn" and
    `> 1.0` as "error" (EFT regime breakdown).

### Changed

- **Three shipped PGT theories re-derived under `[perturbation]`**:
  `torsion_gertsenshtein.json`, `graviton_torsion.json`, and
  `torsion_gertsenshtein_combined.json`. Their loaders now flow
  through the v6 iterative path; `base_spec` demotes the previously-
  4th-order fields `h_4/h_7/h_9` to algebraic constraints at ε=0.
- `tidal simulate --scheme modal` eligibility check now uses
  `spec.base_spec()` when the JSON has corrections, so operators
  that appear only in correction sources (d4_t, mixed_T4_*) do not
  cause auto-selection to reject modal.
- `ExportJSON.wl` emits `order_in_eps` on each OperatorTerm and a
  `perturbation` sub-dict in the JSON metadata when
  `[perturbation]` is configured in the theory.toml.

### Removed

- **`tidal/symbolic/ostrogradsky.py`** (547 lines) — legacy
  auxiliary-field reduction. Closes issues #195, #196.
- **`tidal/wolfram/PerturbativeReduction.wl`** (525 lines) — v5 JLM
  symbolic reducer, superseded by order-tagging.
- **`fields=` kwarg** on `compute_system_energy`,
  `compute_energy_timeseries`, `compute_conversion_probability`,
  `compute_group_conversion`. Only existed to work around the
  Ostrogradsky auxiliary-field layout; no longer needed.

### Breaking

- JSONs derived from theories with `time_order > 2` but without a
  `[perturbation]` section in `theory.toml` now raise on load with
  a clear migration hint. Users must add `[perturbation]` and
  re-derive.

### Fixed

- `solve_modal` with `return_eigendata=True` now exposes per-block
  V / D / V⁻¹ / α plus the Schur recovery operator so the Pass 1
  Duhamel can reuse Pass 0's eigendecomposition (no re-work).
- `_compute_validity` uses `max(|λ|)` instead of `max(|Im(λ)|)` so
  damped/tachyonic eigenvalues register at their full magnitude.

### Deferred

- Euler-Heisenberg example (Phase 6C.2) blocked by xAct
  `Validate::repeated` on (F·F)² — tracked in issue #271.

## [0.31.0] — 2026-04-14

### Fixed

- **#256 — Modal solver unified; rank-deficient mass matrices handled
  correctly**. The modal solver previously had two competing matrix
  builders: a `build_constraint_eliminated_matrices` (constraint path)
  that silently mis-classified `d2_t` RHS cross-couplings as spatial
  identity operators (because it used `_EXACT_MULTIPLIERS[op]`, which
  returns only the spatial part — `ones_like(k) = 1` for `d2_t`), and
  a `_build_generalized_evolution_matrices` (generalized path) that
  correctly classified RHS terms by `time_order` and handled rank-
  deficient mass matrices via per-mode eigendecomposition Schur
  elimination. The router preferred the buggy constraint path for any
  spec with `has_constraints=True`. Result: PGT-style models that
  linearize a higher-rank tensor but only use a low-rank projection in
  the Lagrangian (e.g. `torsion_dark_photon` with the `TorsionCDT`
  trace) produced spurious tachyonic modes with `Re(λ) ≈ 1` and
  diverged at long `t_end`. Unified the two builders into a single
  `_build_evolution_matrices` that handles algebraic constraints,
  rank-deficient mass matrices, and near-singular `I − vel_coupling`
  all in one pass; `_evolve_per_mode` dispatches between
  `np.linalg.eig(A)` and `scipy.linalg.eig(A, B)` (QZ decomposition)
  based on whether the builder produced a `B_lhs` matrix. Pinning
  regression tests in `tests/test_solver_modal_unified.py`. See
  `docs/tex/troubleshooting.tex` §"Rank-Deficient Kinetic from
  Trace-Projection Lagrangians" for the diagnostic pattern and
  resolution.

### Changed

- `solve_modal` router simplified to a single unconditional reduction
  branch (`needs_reduction = has_constraints or has_time_ops`) when
  either algebraic constraints or time-derivative operators are
  present.
- `tidal/measurement/_stability.py::check_conversion_stability` now
  imports `_build_evolution_matrices` (fixes a pre-existing latent
  NameError — the old code referenced `build_constraint_eliminated_matrices`
  without importing it).
- `docs/tex/modal_solver.tex` function reference table, algorithm
  routing table, and implementation paragraph updated to reflect the
  unified builder. Cross-reference added to the new troubleshooting
  section.

### Removed

- `build_constraint_eliminated_matrices` and its export from
  `tidal/solver/__init__.py` — replaced by the unified
  `_build_evolution_matrices`. ~240 lines of dead code deleted.
- Fundamental-vector rewrite of `examples/torsion_dark_photon/theory.toml`
  and `examples/data/torsion_dark_photon.json` reverted (commits
  `cfaa655`, `f3f242a`, `e41b288`); the canonical `TorsionCDT`-trace
  Lagrangian `(1/κ²) R + α I₃ − ¼ ξ Ftorsion² + δₘ F·Ftorsion − ¼ F²`
  is restored and now works natively through the unified path.
  `examples/torsion_dark_photon/sweep_xi.sh`, `sweep_2d.sh`, and
  `sweep_mc.sh` are the canonical sweep scripts again;
  `sweep_mT2.sh` (from the reverted fundamental-vector rewrite) is
  gone.

### Migration

- Downstream code that imported `build_constraint_eliminated_matrices`
  should switch to `_build_evolution_matrices` and destructure a
  6-tuple `(A_rhs, B_lhs, recovery, v_recovery, c_names,
  orig_to_reduced)` instead of the 5-tuple produced by the old
  function. The only internal caller
  (`tidal/measurement/_stability.py`) was updated.

## [0.30.2] — 2026-04-13

### Documentation

- Add pre-arXiv primary references (Gertsenshtein 1962, Boccaletti 1970,
  Raffelt-Stodolsky 1988) to `docs/references.md` with transcribed key
  equations and cross-references to the TIDAL files that implement them.
- Create `docs/tex/gertsenshtein_proca.tex` covering the plasma-Gertsenshtein
  physics motivation, the perturbation-level Proca term justification, and
  the v0.30.1 validation results.
- Update `docs/tex/pipeline.tex` with the new "Canonicalization of Deferred
  Derived Fields" section documenting the v0.30.0 architectural fix.
- Add `examples/gertsenshtein_proca/README.md` as a navigation aid from
  physics intent → theory.toml → sweep scripts → literature.

### Deliverables

- `examples/gertsenshtein_proca/reproduce_figures.sh` — single-command
  figure reproduction script (HPC-friendly: reuses existing sweep outputs
  by default, `--fresh` to force re-run). Produces three supervisor
  hand-off figures: 1D Raffelt-Stodolsky Lorentzian, 2D resonance map,
  mA² family oscillation comparison.
- `examples/gertsenshtein_proca/sweep_mA2_family.sh` — new time-evolution
  sweep at fixed B₀=0.10 for mA² ∈ {0, 0.0025, 0.005, 0.010, 0.050, 0.500}
  showing the transition from coherent Rabi oscillation to off-resonance
  suppression.

## [0.30.1] — 2026-04-13

### Added

- **Regression tests** (`tests/test_gertsenshtein_h5_regression.py`): 13
  new tests covering `gertsenshtein_proca`, `torsion_gertsenshtein_propagating`,
  and `graviton_torsion` JSONs. Total regression suite: 22 tests.
- **1D fast-validation sweep** (`examples/gertsenshtein_proca/sweep_resonance_1d.sh`):
  40-point mA² scan at fixed B₀=0.10 producing a clean Raffelt-Stodolsky
  Lorentzian.

### Changed

- `examples/gertsenshtein_proca/sweep_resonance.sh` grid refined from
  `mA2=0.0:2.0:12` (Δ=0.18, too coarse to resolve the fall-off) to
  `mA2=0.0:2.0:30 × B0=0.01:0.25:10` (300 sims). HWHM observed at B₀=0.25
  matches theory to 9%.

### Validation

- 1D sweep: HWHM = 0.405 vs theory 2·κ·B₀·k = 0.402 (0.7% agreement).
- 2D sweep: per-B₀ HWHM tracks theory within 9% across B₀ ∈ [0.09, 0.25].
- All 1738 tests pass (2 pre-existing chern_simons deselects).

## [0.30.0] — 2026-04-13 — Unified deferred-field canonicalization

### Removed

Three layered workarounds for what was ultimately one pathology:

- `LagrangianFullExpand` Module block in `_wls_linearize_from_lagrangian`
  (added in 9f496ef for #250 h_5 kinetic-term loss)
- Torsion-deferred `ToCanonical` bypass (added in 51a36fd for #255 mixed
  GR+torsion cross-sector coefficient corruption)
- `lagForCanon = If[ValueQ[LagrangianFullExpand], ...]` conditional in
  `_wls_canonical_phase_b` and its dual gauge-fixing / PertLagTerms
  injections
- Adaptive per-term batching loop (`canonBatchSize` = 20 with cost-based
  scaling) in `_wls_linearize_from_lagrangian`

### Root cause identified

Per-term batched `ContractMetric[ToCanonical[Total[batch]], metric]`
fragments cross-term tensor structure: when a graviton kinetic contraction
like `∂h_{ij}∂h^{ij}` spans multiple `Plus`-terms, ToCanonical sees only
trace-like fragments per batch and collapses h_5 = h_{xy} off-diagonal
structure into the constraint set. The batching existed purely as a
performance optimization and was semantically unsafe for any expression
with cross-term kinetics. h_5 was the symptom we caught; a future theory
with spin-2 off-diagonal couplings would hit the same issue.

The torsion cross-sector problem (#255) is a second, independent
pathology: abstract `Ftorsion` and abstract `F` together trigger global
tensor symmetries in `ToCanonical` that halve non-trace h_{ab} coefficients.
Both pathologies share the same root: **ToCanonical on L^(2) with abstract
deferred tensors is ill-behaved**.

### Added

Single unified canonicalization path:

1. Expand deferred derived fields (`F → d(a)`, `Ftorsion → d(t)`) BEFORE
   canonicalization via `_wls_deferred_field_expand`. With F and Ftorsion
   concrete, no abstract tensors remain for global symmetries to act on.
2. Single-pass `ContractMetric[ToCanonical[l2Raw]]` inside
   `TimeConstrained[..., 15.0, $Aborted]`.
3. For very large theories where single-pass exceeds 15s (e.g. #201
   non-minimal R̃[μν]F with 5000+ terms), fall through to
   `l2Raw = Expand[l2Raw]`. Downstream `DecomposeScalarExpression` then
   canonicalizes per-component via `ToBasis + TraceBasisDummy` at much
   smaller scale, preserving correctness.

### Retained

- `a50b11f` F-deferral *upstream* of xPert. This is an orthogonal
  xPert-layer concern (avoiding pre-xPert cancellation that broke photon
  EOMs in #218), not a ToCanonical concern. It stays.

### Impact

- 175 lines deleted from `tidal/cli/_derive.py`, 3 code paths collapsed
  into 1
- All four tracked deferred-field JSONs (gertsenshtein, torsion_dark_photon,
  graviton_torsion, torsion_gertsenshtein_propagating) re-derive with
  **bitwise-identical coefficients** to v0.29.4 — only `derivation_hash`
  changes
- 1725 tests pass
- End-to-end physics validated: plasma-Gertsenshtein reproduces the
  Raffelt-Stodolsky Lorentzian resonance (validated in v0.30.1)

### Verification trail

- Probe 0.1 (force `has_torsion_deferred=True` for gertsenshtein): 6/6
  h_5 regression tests pass — confirms the bypass path handles base
  gertsenshtein correctly
- Probe 0.2 (single-pass `ContractMetric[ToCanonical[l2Raw]]` after
  deferred expansion): 9/9 regression tests pass for both gertsenshtein
  and torsion_dark_photon — confirms single-pass is sufficient
- Commit: `4fef016` `refactor: unify deferred-field canonicalization, remove #218/#250/#255 workarounds`

## [0.29.4] — 2026-04-13 — Minimal correctness fix

### Fixed

- `LagrangianFullExpand` (from 9f496ef) silently dropped Ftorsion
  expansion for torsion theories, causing propagating torsion modes to
  be lost from the canonical pipeline for the propagating-torsion model.
  Restricted the `LagrangianFullExpand` Module block to non-torsion
  theories so torsion-deferred cases continue to follow the #255 bypass
  path unchanged.
- Re-derived all tracked JSONs against v0.29.4 code so committed JSONs
  are reproducible from the current pipeline (not from an earlier
  codebase state).

### Added

- `tests/test_gertsenshtein_h5_regression.py` — 9 regression tests
  protecting the #218/#250/#255 behaviors. This is the first codified
  trail of what each workaround was defending against.

### Note

This was a contained minimal fix. The underlying architectural debt
(three layered workarounds for one pathology) was resolved in v0.30.0.

### Commit: `7da64c7`

## [0.29.3] — 2026-04-12 — #250 h_5 minimal workaround

### Fixed

- Issue #250: h_5 = h_{xy} graviton off-diagonal mode was being
  mis-classified as a constraint (time_order=0) instead of a dynamical
  wave field (time_order=2). Root cause: per-term batched `ToCanonical`
  on abstract `F·F` invariants produced a canonical form where graviton
  kinetic terms were trace-only, losing the non-trace `∂h_{ij}∂h^{ij}`
  structure.

### Workaround (superseded in v0.30.0)

Introduced `LagrangianFullExpand` — a second copy of L^(2) with deferred
fields fully expanded BEFORE ToCanonical — processed in parallel with
the main Lagrangian and used only by the canonical pipeline
(`lagForCanon`). This was a local workaround, not a root-cause fix.

### Commit: `9f496ef`
