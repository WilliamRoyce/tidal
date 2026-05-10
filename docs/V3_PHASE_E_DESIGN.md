# Phase E — Localised wavepacket + localised B-field geometry (DEFERRED)

**Created:** 2026-05-10
**Status:** DEFERRED — gated on Phase B (plane-wave architecture) convergence
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)

## Problem

Phase A/B ship with the v2 geometry: plane-wave IC at the source field + spatially-uniform B₀ background. In tachyonic parameter regions, the linearised PDE has eigenvalues with positive real part, so P_max grows as `~exp(γ·t_end)` without bound. Two consequences:

1. The chain's P_max ranking in tachyonic regions is dominated by exp(γ·t_end), making sup chains hard to interpret as physical "minimum conversion".
2. Larger t_end gives larger P_max for tachyonic samples, so the headline A_max is t_end-conventional rather than physical.

The Phase 6 boundary-attractor analysis (`docs/PHASE_6_COMPARISON.md` §"Phase 6.L") quantified this for the v2-gated chains: A_chain(t_end) = A_static × exp(γ_eff·t_end) with γ_eff ≈ 0.137 at the v2 D1 amp MAP, giving an exp(1.37) ≈ 4× residual growth contribution on top of the t-independent A_static ≈ 10.

## Physical resolution: localised geometry

Replace the plane-wave + uniform-B₀ setup with:

- **Localised wavepacket IC**: `--ic gaussian --ic-component h_5 --ic-amplitude 1e-2 --ic-width σ_w --ic-center x_c`. Already supported in the codebase.
- **Localised B-field profile**: convert the constant `B0` parameter into a position-dependent background field via `[[background_fields]]` in the relevant `theory.toml` files. Gaussian profile in x: `B0 * exp(-((x - x_B) / σ_B)**2)`. Width σ_B becomes a new parameter.

Result: the wavepacket traverses the B-field interaction region in finite time `t_int ~ (σ_w + σ_B) / c`. After traversal, source/target amplitudes propagate in vacuum where the EOM is the free wave equation (no Gertsenshtein conversion). Tachyonic eigenvalues no longer accumulate growth indefinitely — total amplification is bounded by exp(γ·t_int), which is physically meaningful (a finite interaction time turns the unbounded exponential into a finite multiplier).

## Tasks

### E.1 — Wavepacket IC switch

Switch all v3 campaign scripts in `scripts/hpc_submit_drafts/v3_permissive/` (and `v3_jointprior/` once it lands) from `--ic plane-wave` to `--ic gaussian`. Add explicit `--ic-width σ_w` and `--ic-center x_c`. Reasonable defaults at L=100 box: `σ_w = 5` (so the wavepacket occupies ~10% of the domain), `x_c = 25` (a quarter into the box, room to traverse before periodic wrapping).

Already supported per CLAUDE.md — no code changes needed.

### E.2 — Localised B-field background

Modify the `[[background_fields]]` block of the relevant theory TOML files (`examples/torsion_gertsenshtein/`, `examples/torsion_gertsenshtein_nonminimal/`, etc.) to declare `B0` as a position-dependent background:

```toml
[[background_fields]]
name = "B0"
profile = "gaussian"
amplitude = 0.01
center = 75.0
width = 25.0
```

(Profile syntax may need extension if `gaussian` isn't already a recognised position-profile type; check `tidal/symbolic/_derive.py` and `tidal/solver/coefficients.py`.)

Re-derive the affected theories via `tidal derive` (Wolfram pipeline). Expected wall: ~30 min per theory.

Width σ_B = 25 (a quarter of the box) at L=100 keeps the wavepacket-B-field overlap region centred. The wavepacket transit time from σ_w = 5 / σ_B = 25 / box L = 100: full traversal in ~50 time units; the existing t_end = 10 is well within "wavepacket inside B-field" regime, so the localised geometry effectively sets a new default t_end based on σ_B / c ~ 25.

### E.3 — Tuning σ_w and σ_B

Brief 2D scan over (σ_w, σ_B) at the D1 v2 amp MAP (`α₁=−0.422, α₂=−0.594, α₃=0.204, δ₁=−0.847`):

- σ_w ∈ {2, 5, 10}
- σ_B ∈ {10, 25, 50}
- Measure A = P_max / P_GR at each (σ_w, σ_B); also at multiple t_end values (5, 10, 20).

Acceptance: A(σ_w, σ_B, t_end) is approximately t_end-independent (plateaus once t_end > σ_w + σ_B), confirming finite interaction time bounds growth.

Output: `examples/data/v3_localised_geometry_tuning/` with the 9-point CSV and a 2D plot.

### E.4 — Re-baseline campaign scripts

Write `scripts/hpc_submit_drafts/v3e_localised/` mirroring `v3_permissive/` but with localised geometry. All 12 paired campaigns runnable on the new geometry.

### E.5 — Re-run highest-value chains under E geometry

Initial scope (smaller than full v3 fan-out):

- D1 amp + sup (publication targets, paired)
- One D2 sub-model amp (D2.0 Bahamonde representative)
- Stage A sup (the structured null result; may shift)

If those land cleanly, expand to D2.1–D2.3.

### E.6 — Compare Phase B (plane-wave) vs Phase E (localised)

Comparison table in `docs/PHASE_E_GEOMETRY.md` (new). Per-coupling marginal D_KL shifts; per-coupling posterior-shape comparison; A_max headline number under each geometry. Decision: which geometry gets the publication numbers.

If shifts are < 0.5 nats marginal D_KL per coupling: Phase B numbers are defensible without the geometry change (publication uses plane-wave with a footnote). If shifts are larger: Phase E geometry becomes canonical, Phase B numbers become an appendix.

## Done condition

`docs/PHASE_E_GEOMETRY.md` written; manuscript decision made on which geometry features in headline numbers; remaining v3 campaign chains either re-run under E or explicitly footnoted as plane-wave.

## Why deferred

Phase B chains under plane-wave geometry establish the v3 baseline numbers and validate the new architecture (soft penalties, compactified priors, no Hwang-Noh). Until those land, there's no baseline to compare Phase E against — the geometry pivot's *value* is precisely the shift it introduces. Running E first would lose that comparison.

Also, Phase E introduces new tuning parameters (σ_w, σ_B) that should be physically motivated. The tuning step (E.3) is a small project in itself, sensible to do once the rest of the architecture is stable.

## References

- [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md) §"Why v3?" #6
- [docs/PHASE_6_COMPARISON.md](PHASE_6_COMPARISON.md) §"Phase 6.L" — boundary-attractor finding that motivates this phase
- [docs/meetings/2026-05-08_supervisor.md](meetings/2026-05-08_supervisor.md) §3
- CLAUDE.md "Background fields" + "Common pitfalls" — `--ic gaussian` and `[[background_fields]]` are already supported
