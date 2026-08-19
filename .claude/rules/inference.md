---
paths:
  - "tidal/inference/**"
  - "scripts/analysis/**"
---

# Inference Architecture Rules

## Marginal D_KL invariants (#420 family — violating any of these reintroduces the bug)

1. **Every marginal is estimated in the space where its prior is uniform.**
   `_uniformizing_transform` in `tidal/inference/_importance.py` must cover EVERY
   kind in `tidal.inference._prior.VALID_DISTRIBUTIONS` — enforced by
   `test_every_valid_kind_has_uniformizing_transform`. Adding a prior kind means
   adding its transform (or an empirical-reference branch) in the same change.
2. **arctan_uniform's recorded low/high are UNUSED** (#425): the support is fixed at
   ±tan(π/2 − `_ARCTAN_EPS`) ≈ ±19.98. Never derive a histogram range, support, or
   density from those recorded bounds. `0:0` is the sanctioned "unused" sentinel.
3. **Never rank a marginal at or below its noise floor** `(n_bins − 1)/(2·n_eff)`
   (#433). Every consumer of `marginal_d_kl` (tables, plots, TeX emitters) must
   consult `consistency["floor_dominated_params"]` / `noise_floor` before ranking
   or quoting. Floor values are estimator bias, not constraint.
4. **The consistency block is load-bearing, not decoration.** importance.json
   `schema_version: 2` carries: superadditivity (only meaningful when
   `superadditivity_applicable`), `n_eff`, `noise_floor`, `fallback_params`
   (scored against the sample range, not a prior — unreliable), `range_clipped`
   (posterior mass outside the recorded prior range — the prior record does not
   match the samples). Bump `schema_version` on ANY estimator or block change;
   a pre-v0.48.8 importance.json has no consistency block and must not be quoted.
5. **Prior provenance is content, not file existence.** An `inference.json` without
   a `priors` key means fabricated priors (see `recompute_parameter_kl.py`
   `_read_priors` → `priors_provenance`); the Barker-amp chain (#434) was scored on
   fabricated arctan priors for months because existence was mistaken for
   provenance. Positional prior/column alignment mismatches must RAISE, never
   proceed with misassigned records.
6. **Silent failure is the disease.** `InferenceResult.save` logs importance
   failures at WARNING; the estimator raises on degenerate weights; fallbacks are
   recorded machine-readably AND logged. If you add a fallback, record it in the
   consistency block.

## Results record policy

- `manuscript/` is a FROZEN archive of the thesis — never edit anything under it
  (including generated listings; `scripts/analysis/extract_d_kl.py` writes there —
  do not run it without explicit user instruction).
- `CAMPAIGN.md` is append-only: corrections are new dated sections, never edits.
- The living corrected-results record is `docs/RESULTS_AMENDMENTS.md`; the #420
  evidence base is `docs/dkl_recompute_report.md`.
