# Modal A-template Cache — Failure Retrospective

**Reverted:** 2026-05-03
**Reverted commits:**
- `62352c4` — perf(modal): pre-Schur ingredient template cache for constraint specs
- `a068453` — chore: bump version to 0.38.5

**Revert commits:**
- `b4e70c6` — Revert "perf(modal): pre-Schur ingredient template cache for constraint specs"
- `f994329` — Revert "chore: bump version to 0.38.5"

To recover the exact broken implementation for forensic comparison:
`git show 62352c4 -- tidal/solver/modal.py tidal/inference/_likelihood.py tests/test_modal_template_cache.py`.

---

## 1. What this was supposed to do

D2.x amplification chains with `nlive=800` need ≤ 15 ms per likelihood call to
fit in a 15-min CSD3 INTR slot. The pre-revert per-likelihood wall on the
`torsion_gertsenshtein_general_nonminimal.json` fixture (34 fields, 745 RHS
terms, 11 coupling parameters) was ~72 ms, dominated by:

| Component | ms/call |
|---|---|
| Probe (all-k unit-IC, eigenvalue check) | ~13 |
| `CoefficientEvaluator.__init__` | ~19 |
| `_build_evolution_matrices` | ~30 |
| `_ifft_slots` / irfft | ~9 |
| `_evolve_per_mode_pade` | ~0.3 |
| **Total** | **~72** |

The plan was to template the four pre-Schur ingredient matrices (A_dd, S_cc,
S_cd, A_dc_*) — each claimed to be exactly linear in free coupling parameters —
build them once at `base_params`, then assemble per-call via `(N_free)` linear
combinations and a fast batch Schur elimination (~3 ms). Target: ~26 ms/call,
~2.8× speedup.

Reality: it was wrong on correctness AND slower in practice. Both directions
matter.

---

## 2. Failure modes

### 2a. The linearity assumption is false

**Claim:** Pre-Schur ingredients are exactly linear in free coupling params
because every coefficient in the Lagrangian is linear in those params.

**Reality:** The four ingredient matrices captured AFTER the M-matrix
eigendecomposition / basis rotation are NOT linear in any param that affects M.
The basis rotation `Uc_T @ K_mat @ V_full` ([tidal/solver/modal.py:1175-1245](
../tidal/solver/modal.py#L1175-L1245)) makes V_full param-dependent whenever
the mass matrix M has a non-trivial dependence on a coupling — which is
precisely the case for `delta1` in the GENERAL spec (kinetic mixing /
non-minimal coupling).

**Empirical evidence** (template built at delta1=0, applied at varying delta1,
GENERAL spec, grid=32):

| param tested | A_dd error at d=0.5 | other ingredients |
|---|---|---|
| delta1 | 0.80 abs (31.6% rel) | S_cc/S_cd/A_dc_field/A_dc_vel: 4e-16 ✓ |
| beta2  | 4.4e-16 (machine ε) | all linear ✓ |
| zeta1  | 0 | all linear ✓ |

**Error scaling** confirms quadratic non-linearity in delta1 (linear-interp
error `a·d·(d−h)` with probe at `h=1`):

| delta1 | A_rhs abs err | A_rhs rel err |
|---|---|---|
| 0.001 | 3.2e-3 | 0.015% |
| 0.010 | 3.1e-2 | 0.15% |
| 0.025 | 7.5e-2 | **0.36%** |

At campaign delta1 ∈ [-0.025, 0.025], the post-Schur A_rhs is wrong by ~0.36%
relative — enough to bias a Bayesian posterior tail. S_cc / S_cd / A_dc_field /
A_dc_vel ARE exactly linear; the corruption is contained to A_dd, but A_dd
contributes directly to A_rhs so it propagates.

### 2b. Spec re-loading defeats the cache

`run_inference_step` ([tidal/cli/_sweep.py:619-658](../tidal/cli/_sweep.py#L619-L658))
calls `load_equation_system(spec_path)` on every likelihood call when
`spec=None` is passed (which is what `SimulationLikelihood.__call__` does).
The cache key in the reverted code was `(id(spec), grid_shape, grid_bounds)` —
since `id(spec)` differs every call, **the cache miss rate is 100%**.

Each "miss" rebuilds the template at cost `(N_free + 1) × full_build_cost` —
i.e. 2× the legacy build cost for `N_free = 1`, scaling linearly with
`N_free`. Then `_apply_ingredient_template` is invoked with
`parameters == base_params` (the freshly-built base), so `delta = 0` and it
just returns A_dd_base unchanged. **Net effect: always pay the build penalty,
never amortize.**

This also explains why the simulation-equivalence sanity test reported
bit-exact agreement: every call rebuilt; `_apply` ran at delta=0; A_dd_base
happened to equal the legacy direct build at the same params.

### 2c. The benchmark proved it was slower

Microbenchmark, GENERAL spec, 10 calls via `run_inference_step`,
delta1 randomized over [-0.025, 0.025]:

| mode | median ms | mean ms | min ms |
|---|---|---|---|
| Template OFF | 82.6 | 83.9 | 64.3 |
| Template ON  | 104.5 | 108.3 | 84.8 |
| **Speedup** | **0.79× (21% SLOWDOWN)** | | |

### 2d. No build-time validation

Nothing in `_build_a_template_cache` verified the linearity assumption at a
midpoint of the prior range. The template silently returned wrong A_rhs.

### 2e. Module-level state pollutes pytest

`_free_param_names`, `_A_TEMPLATE_CACHE`, `_MULTIPLIER_CACHE` were module-level
dicts mutated by `register_inference_context`. State leaked across pytest
sessions — the broader suite tripped
`tests/test_perturbative_driver.py::TestPerturbativeSolverValidity::test_validity_formula_is_eps_times_omega_times_t`
(passes in isolation, fails in full suite). After the revert, this test
ordering issue is gone.

### 2f. Tests that didn't actually test the cache

`tests/test_modal_template_cache.py` had:
- API errors (`EquationSystem.from_file` doesn't exist; correct API is
  `load_equation_system`).
- A stale variable name (`_FREE_PARAM_NAMES`, renamed mid-session to
  `_free_param_names`).
- All 5 tests rebuilt the template fresh per call — they would have passed
  even on the broken implementation, because `_apply` at `delta=0` returns
  A_dd_base which equals the direct rebuild.

The tests were never run during development (they only got
`pytest --collect-only`'d). When fixed and run, the linearity test failed
with 7.7% of A_dd entries differing by up to 0.04.

---

## 3. What was assumed correctly vs. wrongly

### Correct
- The raw coefficient structure is linear in coupling params (no `param_i ×
  param_j` products in the Lagrangian).
- `S_cc`, `S_cd`, `A_dc_field`, `A_dc_vel` are individually linear in coupling
  params (verified to machine precision).
- The Schur elimination is the dominant non-linearity for the *post-Schur*
  output (S_cc_inv is rational in params).
- Stage 1 sparse-IC mode skip ([8ecaf76](https://github.com/WilliamRoyce/torsion-gertsenshtein/commit/8ecaf76))
  is unaffected and works as intended.

### Wrong
- **The original "n_c=0, post-Schur template is exact" claim was false** — but
  it was caught and revised before commit (see plan archive in
  `git show 62352c4` commit body, which already documented this revision).
- **The revised "pre-Schur ingredients are exactly linear" claim was also
  false** — overlooked the M-matrix eigendecomposition. This is what the
  reverted code shipped.
- The performance target of 26 ms/call was measured in the wrong scope. A
  proper end-to-end benchmark via `SimulationLikelihood.__call__` (with the
  spec-reload path) would have shown 2× the legacy build cost up-front.

---

## 4. Prerequisites for any future attempt

Any future template-cache attempt must satisfy ALL FOUR before being merged:

1. **Build-time linearity check at a midpoint of the prior range.** After
   building the template, evaluate at e.g. `delta = 0.5 × (prior_max −
   prior_min)`, compare the post-Schur output to a direct rebuild. If
   `np.allclose(A_rhs_template, A_rhs_direct, atol=1e-10) == False`, return
   None and fall back to legacy. The reverted code had no such check. This
   alone would have caught the failure on day one.

2. **Stable cache key.** Either:
   - Cache the loaded `EquationSystem` in `SimulationLikelihood.__init__`
     and pass it through `run_inference_step` (instead of letting it reload),
     OR
   - Use a content hash of the spec (`hashlib.sha256` over a canonical JSON
     serialization) as the cache key instead of `id(spec)`.

   Without this, the optimization is dead on arrival regardless of correctness.

3. **Microbenchmark with the actual inference call path.** Measure
   `SimulationLikelihood.__call__` end-to-end (or at minimum,
   `run_inference_step`) over ≥ 20 calls with realistic random parameters.
   The 72 ms / 26 ms numbers in the original plan came from a hand-rolled
   `_build_evolution_matrices` benchmark that didn't include spec reload,
   coefficient evaluation, or the cache-miss penalty.

4. **Tests that exercise cache reuse.** Build template once at `base_params`,
   call N times at varying params, assert post-Schur outputs match direct
   rebuild to `atol=1e-10`. The reverted tests rebuilt fresh per call and
   would have passed even on the broken implementation. A test that catches
   this failure mode is non-negotiable.

---

## 5. Reverted artifacts

```
$ git show --stat 62352c4
 tests/test_modal_template_cache.py | 453 +++++++++++++++++++++++++++++++
 tidal/inference/_likelihood.py     |   3 +
 tidal/solver/modal.py              | 540 +++++++++++++++++++++++++++++++++----
 3 files changed, 945 insertions(+), 51 deletions(-)
```

**Symbols removed from `tidal/solver/modal.py`:**
- `_ATemplateCache` dataclass
- `_IngredientTemplateCache` dataclass
- `register_inference_context()`
- `_make_template_cache_key()`
- `_apply_a_template()`
- `_apply_ingredient_template()`
- `_build_a_template_cache()`
- Module-level globals `_free_param_names`, `_A_TEMPLATE_CACHE`,
  `_MULTIPLIER_CACHE`
- Optional kwarg `_ingredient_out` on `_build_evolution_matrices`
- Environment kill-switch `TIDAL_MODAL_A_TEMPLATE`

**Removed from `tidal/inference/_likelihood.py`:**
- Call to `register_inference_context(param_names)` in
  `SimulationLikelihood.__init__`.

**Removed entirely:**
- `tests/test_modal_template_cache.py` (5 tests, all broken).

---

## 6. What lives on (unchanged by this revert)

- **Stage 1 sparse-IC mode skip** ([8ecaf76](https://github.com/WilliamRoyce/torsion-gertsenshtein/commit/8ecaf76))
  — independent and working. Reduces per-call wall when the IC is sparse
  (single-mode plane-wave). Tests in `tests/test_modal_sparse_ic.py`.
- **Stage 3 multiplier cache** was bundled into 62352c4 and reverted with it.
  It showed no measurable benefit in the benchmark (template-OFF baseline of
  82.6 ms vs the original 72 ms plan estimate suggests the cache wasn't
  actually amortizing anything). It can be revisited independently under
  prerequisite #3 above.

---

## 7. Path forward (out of scope for this revert)

If the inference performance target (≤ 15 ms per likelihood at `nlive=800`) is
still required, the right sequence is:

1. **Fix spec lifecycle first** — make `SimulationLikelihood` cache its
   `EquationSystem` and pass it into `run_inference_step`. This is a
   prerequisite for any matrix cache, but is also independently useful: the
   current code reloads + parses the JSON every call, which is itself ~5–10 ms.
2. **Then re-benchmark** without any matrix cache to see what the actual
   bottleneck is.
3. **Only then** consider a template cache, with all four prerequisites
   above satisfied.

This retrospective and `feedback_modal_template_cache_failure.md` (in the
Claude memory) are the durable artifacts — re-attempting the optimization
naively without addressing the prerequisites will reproduce this same failure.
