# Changelog

All notable changes to TIDAL are documented in this file. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
uses semantic versioning as interpreted by `scripts/bump_version.py`.

This changelog was introduced in v0.30.2 and backfilled for the
canonicalization-architecture trail (v0.29.3 → v0.30.1) because those
commits reworked a load-bearing part of the Wolfram pipeline that future
maintainers need to be able to trace. Earlier versions are not
retroactively covered; see `git log` for the full history.

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
