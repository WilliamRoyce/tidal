# v3.2 design investigation — symmetric priors + PSALTer integration

**Created:** 2026-05-11
**Status:** INVESTIGATION DOC — implementation gated on supervisor classification of kinetic vs mass-like parameters
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md), [V3_PHASE_C_REFERENCE.md](V3_PHASE_C_REFERENCE.md)

## Why this document exists

After Phase B.1 (D1 amp) and B.4a (Stage A amp) v3 chains landed, user review of the corner plots surfaced three coupled concerns that go deeper than tactical fixes:

1. **Corner plot scale cropping** — appeared as hard vertical cuts on cross-panels
2. **Asymmetric priors** — kinetic-like couplings (α₃, mA², ξ) use `log_uniform:1e-3:1e3` (positive-only); sign-symmetric couplings (α₁, α₂, δ₁) use `arctan_uniform:-89:89`. The v3 thesis is "remove physical assumptions, let the data speak" — but we're still imposing kinetic positivity by prior fiat. This BLOCKS uniform application of the cubed-sphere joint prior across all couplings
3. **Ghost vs tachyon distinction** — letting kinetic couplings go negative may produce physically meaningless results (Hamiltonian unbounded below); detecting ghosts properly is a hard problem that the supervisor's PSALTer package already solves

This document captures the investigation outcomes and decisions so v3.2 implementation can proceed cleanly once gating questions are resolved.

## Finding 1 — Corner plot cropping is mostly honest

On Stage A's `deltam`-vs-`logA` cross-panel, contours showed a sharp vertical edge at `deltam ≈ +8`. Direct inspection of the weighted chain CSV shows **99.8% of posterior weight is at `deltam ≤ +8.5`**, with chain samples extending to ±20 but carrying negligible weight in the tail.

This is anesthetic's 95%-credible contour correctly stopping at the iso-density boundary — a *real* posterior cliff, not a visualization defect. The genuine physics is a sharp logL falloff above some deltam threshold (likely a stability boundary in the plasma model).

**Improvement opportunity (cosmetic)**: add outer 99%/99.9% credibility contours so the tail-falloff is visible alongside the 95% contour. Single-file change in `tidal/inference/_visualize.py` plot_2d call (`levels=` argument). Optional, not blocking.

Filed as GH issue `[v3-viz] Add outer 99% contour to corner plots for tail visibility`.

## Finding 2 — Joint-prior + positive-only-kinetic: four candidate approaches

| Approach | Math | LOC | Trade-offs |
| --- | --- | --- | --- |
| **A: Per-coupling positivity fold** | After `c = r·θ̂`, apply `c_i = r·\|θ̂_i\|` *only on dims marked positive-only* in the prior spec; sign-symmetric dims keep `c_i = r·θ̂_i` unchanged | ~15 | Same `r ~ log_uniform` magnitude scaling across all couplings; sign-symmetric dims keep full sphere direction freedom; only kinetic dims (whichever are correctly classified) get their negative-θ̂ reflected to positive; quick; sample-time enforcement |
| **B: Log-space reparameterisation** | Sphere over `(log\|α_kin\|, α_other, ...)`; `c_i = ±exp(r·θ̂_i)` for log-space dims | ~25 | Sphere stays uniform on log-coords; BUT magnitude `r` rescales exponentially for log-space dims while sign-symmetric dims see `r` linearly — **inconsistent magnitude treatment across couplings**; hidden exponential transform |
| **C: Split into separate joint priors** | Two `RadialAngularPrior` instances: one sign-symmetric, one positive-only | ~40 | Fragments the supervisor's "monolithic sphere" directive |
| **D: Full symmetry + PSALTer-tagged ghosts** | Sample all couplings sign-symmetrically; PSALTer post-hoc tags `run_status="ghost"` | ~280 (PSALTer bridge) | Mathematically purest — chain genuinely explores all signs; D_KL reveals positivity as a *learned* constraint; requires PSALTer integration to land first |

**Decision**: **Approach A is the v3.2 interim path; Approach D is the ultimate goal**.

- Approach B is rejected because the "hidden log" gives positive-only couplings a different magnitude scaling than the sign-symmetric ones, defeating the purpose of `r ~ log_uniform` as the single scale for all parameters.
- Approach C is rejected because it fragments the supervisor's "monolithic sphere over all BSM couplings" directive (per `docs/V3_PHASE_C_REFERENCE.md` line 61).
- Approach D requires PSALTer integration which is non-trivial; defer to Phase C+ / Phase E+.
- Approach A is minimal (~15 LOC) and preserves the monolithic-sphere geometry. Sign-symmetric dims keep full sphere direction freedom; only kinetic dims get folded.

## Finding 2b — Which parameters are *actually* kinetic? (user-flagged gap)

The v2 priors marked these as `log_uniform` (positive-only):

- **α₃** in D1 (torsion_gertsenshtein_nonminimal)
- **mA², ξ, α₃** in Stage A (dark_photon_plasma)

I assumed this carried over to v3 because they were "kinetic-like". **This is a v2-inherited assumption that needs theory-by-theory justification, not a blanket rule.**

The distinction that matters in v3:

- **Mass-like term** (e.g. `m² φ²` in the Lagrangian): negative coefficient → **tachyon** → v3 *allows* (the whole point of the v3 pivot)
- **Kinetic-like term** (e.g. `(∂φ)²` in the Lagrangian): negative coefficient → **ghost** → must be excluded (unbounded-below Hamiltonian)

Concrete examples:

- **mA² in Stage A**: looks like a mass-squared term `mA²·A_μA^μ`. Negative mA² is a tachyon, not a ghost. **Should be sign-symmetric in v3, not log_uniform-positive.** Currently the v2-inherited prior is wrong.
- **α₃ in D1 vs α₃ in Stage A**: same symbol, possibly different operator. In D1 (torsion-Gertsenshtein non-minimal) α₃ is one of the dimensionless torsion-quadratic invariant coefficients. In Stage A (dark-photon-plasma) α₃ is a torsion-photon coupling. Different operators, different ghost/tachyon implications.
- **ξ in Stage A**: needs identification — kinetic coefficient, mass-mixing parameter, or dimensionless coupling?

**This affects observed B.4a/B.4b results**: Stage A sup (29199129) MAP has `alpha3 = 0.00117`, hitting the v3 prior lower bound `1e-3`. The chain wants to push α₃ even smaller (or potentially negative if it's mass-like), being clipped by prior. The strong joint D_KL = 7.30 nats may be partly driven by this prior-edge effect.

**Action: supervisor clarification needed per parameter per theory** before we can correctly apply Approach A. Filed as GH issue `[v3.2-classify] Kinetic vs mass-like per-parameter classification (supervisor input)`.

## Finding 3 — PSALTer ghost-detection integration

PSALTer (supervisor's package, arXiv 2406.09500 v1 and 2506.02111 v2) exposes:

- `ParticleSpectrum[Lagrangian, options]` (Wolfram/xAct frontend) — saturated propagator + pole-residue analysis
- `psalter._classify.core.classify()` (Python-native, JAX-compiled, microseconds per sample after one-time compilation)
- v2 algorithmic advance: ghost detection via kinetic-matrix definiteness, bypassing radicals for parity-violating theories

Bridge requires:

- TIDAL JSON → PSALTer theory object converter (~200 LOC, new module)
- Evaluator instantiation at likelihood initialization (~50 LOC)
- Per-sample residue inspection + `run_status="ghost"` tagging (~30 LOC)

Why this path over a homegrown check: PSALTer correctly handles gauge-fixing + constraint structure + parity-violating cases (chequer-Hermitian blocks). Naive bare-kinetic-coefficient sign checks are unreliable for PGT theories.

**Target timeline: Phase C+ / Phase E+, not now**. The integration is a 1–2 day implementation IF PSALTer's loader is stable, but waiting for our JSON format and Phase C joint-prior conventions to land first is prudent.

Filed as GH issue `[v3.2-psalter] PSALTer integration for chain-sample-time ghost tagging`.

## Decisions

1. **Don't cancel current Phase B chains** (B.1, B.2, B.3, B.4a, B.4b). They stand as the v3.1 baseline with documented asymmetric per-param priors. They have genuine scientific value and validate the v3 methodology shift.
2. **No homegrown ghost detection in TIDAL**. User explicitly defers ghost-detection to PSALTer integration.
3. **No prior architecture changes until investigation complete**. Continue Phase B as-is under v3.1 priors; v3.2 is a follow-up that happens only after supervisor reviews v3.1 results AND classifies which params are kinetic vs mass-like.
4. **Approach A (per-coupling abs-fold) is the v3.2 interim plan**; Approach D (PSALTer-tagged) is the ultimate goal. v3.2 implementation is GATED on `[v3.2-classify]` issue completion.

## Supervisor question (draft for next meeting)

> The v3 architecture allows tachyonic samples (negative mass-squared) but excludes ghosts (negative kinetic). We currently treat all our v2-inherited `log_uniform` priors as positive-only — α₃ in D1, plus mA², ξ, α₃ in Stage A. But mA² looks like a mass-squared term that *should* be sign-symmetric under the v3 tachyon-permissive policy (negative mA² = dark-photon tachyon, allowed).
>
> Could you confirm per parameter per theory: which are kinetic-coefficient (positivity required, ghost risk) vs mass-like (sign-symmetric OK in v3)? Once we know this, we can apply the per-coupling abs-fold correctly in the cubed-sphere joint prior, OR wait for PSALTer integration which would determine this automatically from the kinetic-matrix definiteness.
>
> Concretely, the Stage A sup chain just landed with alpha3 MAP ≈ 0.0012 hitting our log_uniform lower bound of 1e-3 — the chain wants to go smaller. If alpha3 is actually mass-like in this theory, sign-symmetric sampling would let it cross zero and we'd discover the genuine constraint from the data.

## Out of scope (explicit)

- **Implementing** Approach A (per-coupling abs-fold) — gated on `[v3.2-classify]` supervisor review
- **Implementing** PSALTer integration (Approach D) — gated on Phase C convergence
- **Re-running** v3.1 chains under symmetric priors — gated on v3.2 design completion
- **Homegrown ghost detection** in TIDAL — explicitly user-deferred to PSALTer

## References

- [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md) — canonical v3 architecture
- [V3_PHASE_C_REFERENCE.md](V3_PHASE_C_REFERENCE.md) — cubed-sphere joint-prior implementation
- [literature/2406.09500/](../literature/2406.09500/) — PSALTer v1 paper
- [literature/2506.02111/](../literature/2506.02111/) — PSALTer v2 paper (parity-violating extension)
- `psalter.tar.gz` (gitignored, read via `tar -xzOf` only) — supervisor's unpublished Python evaluator
