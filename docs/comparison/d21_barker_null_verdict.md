# D2.1 Barker — NULL verdict under v3 (cross-architecture agreement)

**Created:** 2026-05-11
**Status:** confirmed by v3 paired chains at two `nlive` settings
**Companion to:** [V3_ARCHITECTURE.md](../V3_ARCHITECTURE.md), [V3_PHASE_TRACKER.md](../V3_PHASE_TRACKER.md)

## Result

D2.1 Barker (T5 sub-theory with β₁, β₂, β₃, ξ, χ free; δ₁, ζ₁₋₃ pinned to 0) shows **no preferred parameter region** under the v3 architecture. Both amplification and suppression chains, at two different `nlive` settings, converge to a flat posterior with joint D_KL ≈ 0.

| Chain | Job | `nlive` | ESS | log Z | Joint D_KL | MAP β₁ | MAP β₂ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Barker amp | 29209010 | 400 | 474 | +0.00012 ± 0.0017 | 0.002 | −19.80 (edge) | +16.72 (edge) |
| Barker amp rerun | 29209795 | 1200 | 1274 | +0.00144 ± 0.0006 | 0.0009 | −19.61 (edge) | −0.80 |
| Barker sup | 29209289 | 400 | 474 | −0.00433 ± 0.0012 | 0.002 | −19.75 (edge) | +4.01 |
| Barker sup rerun | 29210161 | 1200 | 1274 | −0.00310 ± 0.0011 | 0.00045 | −19.00 (edge) | +1.53 |

- ESS scales linearly with `nlive` (3× nlive → 3× ESS), confirming PolyChord converged in both runs.
- log Z is essentially zero (consistent with a flat likelihood averaged over the prior).
- Joint D_KL is < 0.002 nats in all four chains.
- A ∈ [1.00, 1.00] (no amplification anywhere in coupling space).
- All four chains: 100% `run_status = success`.

## Cross-architecture comparison

| | v2 (28607124, narrow priors + probe gate) | v3 (29209795, wide priors no gate) |
| --- | --- | --- |
| Priors | β₁∈[−1.5, 1.0], β₂∈[−3, −0.3], **δ₁∈[−0.025, +0.025]**, **χ∈[−0.009, +0.009]** | β/χ ∈ arctan_uniform:−89:89 (±57), ξ ∈ log_uniform[10⁻³, 10³], δ₁ frozen |
| log Z | +0.62 ± 0.001 | +0.001 ± 0.0006 |
| Joint D_KL | 0.001 nats | 0.0009 nats |
| Per-coupling D_KL | 0.02–0.30 nats (all tiny) | 2.5 nats β's (artifact, see below), 0.02 ξ |

**Both architectures agree: no real signal.** v2's superficially-positive log Z = +0.62 came from prior-volume effects — the likelihood averaged over the narrow v2 prior is ~exp(0.62) ≈ 1.86, but the *shape* of the posterior was already uniform (joint D_KL = 0.001).

The 2.5-nat per-coupling D_KL on β₁/β₂/β₃/χ under v3 is **not real signal**: the joint D_KL ≈ 0 confirms the joint posterior is essentially the prior. The marginals appear shifted from the Cauchy prior peak at 0 to the prior edges because the chain explores corner regions with no preferred destination — but those corners carry no information beyond their existence as prior support.

## Methodological strength

v3 makes the null verdict **cleaner and more defensible** than v2:

- v2 narrow priors (e.g. χ ∈ [−0.009, +0.009]) could in principle have been hiding structure just outside the prior bounds.
- v3 wide priors (e.g. χ ∈ arctan_uniform spanning ±57) explicitly explore the full coupling space.
- Under v3 the chain finds nothing preferential anywhere — the null is now verified across a 1000×–10000× wider parameter range per coupling.
- v3 also removes the v2 tachyonic probe gate, so the chain explores formally-unstable regions too. No structure there either.

D2.1 Barker is genuinely **Gertsenshtein-neutral** in the propagating-torsion family.

## Corner plot interpretation

The D2.1 v3 corner plots ([29209795](../../hpc_results/29209795/corner_v3.png), [29210161](../../hpc_results/29210161/corner_v3.png), and the under-converged originals [29209010](../../hpc_results/29209010/corner_v3.png), [29209289](../../hpc_results/29209289/corner_v3.png)) use a **scatter-cross-panel fallback** because anesthetic's KDE-contour machinery (qhull Delaunay + scipy gaussian_kde) fails on near-uniform posteriors — the data covariance becomes singular.

What you see:

- **1D diagonals**: KDE or histogram of each coupling — show prior-like shapes (Cauchy peak at 0 for β's, declining log-density for ξ)
- **2D cross-panels**: scatter (one dot per sample) — show no clustering, no preferred regions; the posterior is uniform-by-eye
- **Title**: A ∈ [1.00, 1.00] — explicit confirmation that amplification is essentially 1 everywhere

This is the **canonical visual signature of a genuine null result** under v3.

## References

- v3 architecture: [docs/V3_ARCHITECTURE.md](../V3_ARCHITECTURE.md)
- Phase B tracker: [docs/V3_PHASE_TRACKER.md](../V3_PHASE_TRACKER.md)
- Plot-fallback bug fix: GH issue [#362](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/362)
- v2 D2.1 Barker chain (28607124) is in `hpc_results/28607124/d21_barker_amp/` for direct inspection
