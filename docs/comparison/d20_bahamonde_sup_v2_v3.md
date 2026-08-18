# D2.0 Bahamonde sup: v2 vs v3 comparison

**Status:** Truncated result — supervisor-approved 1-INTR-session policy (2026-05-15).
The v3 chain was canceled after 7 INTR sessions (job 29271186 canceled at session 8 boundary,
last completed checkpoint from job 29269357 session 7). Not converged; logZ still changing.

**Corner plot:** `hpc_results/29232780/corner_d20_sup_v3.png`

## Key numbers

| Metric | v2 (28519675) | v3 (29232780, truncated) |
|---|---|---|
| logZ | −0.449 ± 0.002 | +134.47 ± 0.453 (not converged) |
| ESS | 877 | — (chain not converged) |
| Dead points | ~full | 58,457 |
| Active clusters | 0 (converged) | 174 / 401 total |
| Run status | Standard 12 h queue | 7× INTR (canceled) |

## Parameters

| Param | v2 MAP | v3 prior | v3 dead-pt mean |
|---|---|---|---|
| beta1 | −0.226 | arctan_uniform[−89°..89°] | −0.43 |
| beta2 | −1.340 | arctan_uniform[−89°..89°] | −0.61 |
| beta3 | −0.294 | arctan_uniform[−89°..89°] | −1.15 |
| xi | 0.211 | log_uniform[1e-3..1e3] | ~649 |
| delta1 | −0.024 | arctan_uniform[−89°..89°] | −0.28 |

**Note:** v3 dead-pt means are from `ns.compress()` over 58K dead points while chain
is not converged; they should be treated as landscape indicators, not posterior estimates.

## Interpretation

- **logZ shift of +134.92 nats is not a like-for-like comparison.** v3 uses
  arctan_uniform priors (effectively ±57 in coupling units) versus v2's narrow
  uniform priors (typically ±3 or ±0.3). The prior volume expansion alone accounts
  for a large fraction of the shift; the tachyon-permissive architecture then adds
  contributions from previously-rejected parameter regions.
- **401 total clusters** (vs 27 for the amp chain at D2.0 Bahamonde) confirms that
  the suppression likelihood landscape is dramatically more multimodal than amplification.
  Many isolated suppression modes each contribute separately to Z.
- **No convergence signal from this run.** σ_logZ was still 0.453 at truncation
  (decreasing from 0.460 — some convergence onset — but 174/401 active clusters
  still active). Under the new 1-INTR-session landscape-overview policy, this state
  is accepted as the v3 sup result.
- **xi MAP jumped to ~649** (v2: 0.211) — a massive prior-shift driven by
  tachyon-permissive regions at large xi that v2 rejected outright.

## Policy note (2026-05-15 supervisor meeting)

PolyChord chains for v3 campaign are treated as landscape-overview runs, not
convergence-required results. Each chain runs 1 INTR session (≤1 h); state at
TIMEOUT is accepted. This chain consumed 7 sessions and is therefore an outlier;
going forward, D2.1–D2.3 will each receive exactly 1 INTR session.
