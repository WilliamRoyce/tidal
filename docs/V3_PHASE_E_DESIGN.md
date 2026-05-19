# Phase E — Localised wavepacket + localised B-field geometry

**Created:** 2026-05-10
**Status:** E.0 empirical validation complete (2026-05-16). E.1–E.6 still deferred pending Phase B convergence.
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

## E.0 — Dual-Gaussian localised B-field (empirically validated 2026-05-16)

### Why not single-Gaussian

The original sketch in this document (§"Physical resolution") proposed a single-Gaussian B₀ profile. The integral of a single Gaussian is an error function:

```text
A_y(z) = -Bpeak · R · √(π/2) · Erf[z/(√2·R)]
```

which runs from −1 to +1 across the domain. This is explicitly non-periodic: `A_y(0) ≠ A_y(L)`. Consequence: any theory TOML using this form must declare **Neumann BCs** (as `examples/gertsenshtein/theory_localized.toml` does explicitly on line 24), which prevents `can_use_modal()` from returning True — the modal solver is never auto-selected for single-Gaussian localised-B runs.

### Dual-Gaussian trick (supervisor, 2026-05-15)

Use two equal-and-opposite Gaussian peaks with centres symmetric about L/2:

```text
B_z(z) = Bpeak · [exp(-(z-zc1)²/(2·sigB²)) − exp(-(z-zc2)²/(2·sigB²))]
```

with `zc1 + zc2 = L`. Then `∫₀ᴸ B_z dz = 0`, and the gauge potential

```text
A_y(z) = -Bpeak · sigB · √(π/2) · [Erf[(z-zc1)/(√2·sigB)] − Erf[(z-zc2)/(√2·sigB)]]
```

satisfies `A_y(0) = A_y(L) = 0` (both Erf terms saturate to ±1 at the boundaries; the difference cancels). With `sigB ≪ zc1`, boundary leakage `B(0) ~ exp(-zc1²/(2·sigB²))` is negligible — at `sigB=5, zc1=25`: `exp(-12.5) ≈ 4×10⁻⁶·Bpeak`.

### Canonical theory file

`examples/gertsenshtein/theory_e0_dual_gaussian.toml` — new (2026-05-16). Key `[[background_fields]]` block:

```toml
[[background_fields]]
name = "Abar"
type = "vector"
components = [
  "0",
  "0",
  "-Bpeak * sigB * Sqrt[Pi/2] * (Erf[(z[] - zc1)/(Sqrt[2]*sigB)] - Erf[(z[] - zc2)/(Sqrt[2]*sigB)])",
  "0"
]

[constants]
names = ["kappa", "Bpeak", "sigB", "zc1", "zc2"]
```

Same fields (h, A, a, F), same Lagrangian, same TT+Lorenz gauges, same plane-wave reduction as `theory_localized.toml`. No Neumann BC override — the TOML uses default periodic BCs.

### Validated parameters

| Parameter | Value | Rationale |
| --- | --- | --- |
| `sigB` | 5 | Decay to 4e-6 at boundaries; smooth-periodic |
| `zc1` | 25 (= L/4) | Left peak at quarter-box |
| `zc2` | 75 (= 3L/4) | Right peak; zc1+zc2=L=100 |
| `t_end` | ≤ 25 | Wavepacket at zc1=25 reaches zc2=75 at t≈50; keep below L/4=25 to avoid wrap-around |
| `--ic gaussian --ic-width 5 --ic-center 25` | | Wavepacket at left interaction region |
| `--grid-shape` | 256 | Resolves sigB=5 at L=100 |

### Boundedness verification (A(20)/A(10) criterion)

Measured with `--scheme cvode`, grid=256, bounds=0:100, Bpeak=0.1, sigB=5, zc1=25, zc2=75, kappa=1.0:

| t_end | a_1 peak | P_max |
| --- | --- | --- |
| 5 | 0.0001 | ~0.0001 |
| 10 | 0.0003 | ~0.0003 |
| 15 | 0.0003 | ~0.0003 |
| 20 | 0.0003 | 0.00348 |

**A(t=20) / A(t=10) = 1.00 ≪ 1.05 criterion ✓** — bounded interaction time, no tachyonic accumulation.

Full conversion: P_max = 0.003475 at t = 18.0 s (then decays as wavepacket exits interaction region).

### Modal solver caveat (GH #367)

The dual-Gaussian theory has position-dependent coefficients (Erf background). `can_use_modal()` still returns True (periodic BCs + flat metric + supported operators). However, `_has_position_dependent_terms()` returns True, routing dispatch to `_evolve_full_matrix()`, which calls `expm_multiply` on a non-normal convolution matrix. The convolution matrix has positive-real eigenvalues (max ~+0.5 at N=16) — a pseudospectral artifact — causing `expm_multiply` to diverge catastrophically (h_5 → 10⁶⁰ by t=20).

**Workaround: always use `--scheme cvode` for Phase E runs.** CVODE gives correct bounded results. Issue #367 tracks the fix.

## Tasks

### E.1 — Wavepacket IC switch

Switch all v3 campaign scripts in `scripts/hpc_submit_drafts/v3_permissive/` (and `v3_jointprior/` once it lands) from `--ic plane-wave` to `--ic gaussian`. Add explicit `--ic-width σ_w` and `--ic-center x_c`. Reasonable defaults at L=100 box: `σ_w = 5` (so the wavepacket occupies ~10% of the domain), `x_c = 25` (a quarter into the box, room to traverse before periodic wrapping).

Already supported per CLAUDE.md — no code changes needed.

### E.2 — Localised B-field background

Create new theory TOML files for each model (`examples/torsion_gertsenshtein/`, `examples/torsion_gertsenshtein_nonminimal/`, etc.) following the dual-Gaussian pattern validated in E.0. The `[[background_fields]]` block uses the explicit Erf formula — there is no `profile = "gaussian"` shorthand in the TOML schema; the position-dependent expression goes directly into `components`:

```toml
[[background_fields]]
name = "Abar"
type = "vector"
components = [
  "0",
  "0",
  "-Bpeak * sigB * Sqrt[Pi/2] * (Erf[(z[] - zc1)/(Sqrt[2]*sigB)] - Erf[(z[] - zc2)/(Sqrt[2]*sigB)])",
  "0"
]

[constants]
names = [..., "Bpeak", "sigB", "zc1", "zc2"]
```

Re-derive the affected theories via `tidal derive` (Wolfram pipeline). Expected wall: ~30 min per theory.

Use `sigB = 5` (validated in E.0) with `zc1 = L/4`, `zc2 = 3L/4`. Cap `t_end ≤ 25` to keep the wavepacket inside the interaction region before it wraps. All Phase E simulations must use `--scheme cvode` until GH #367 is resolved (modal solver auto-selects but diverges for position-dependent Erf backgrounds).

### E.3 — Tuning σ_w and σ_B

Brief 2D scan over (σ_w, σ_B) at the D1 v2 amp MAP (`α₁=−0.422, α₂=−0.594, α₃=0.204, δ₁=−0.847`):

- σ_w ∈ {2, 5, 10}
- σ_B ∈ {10, 25, 50}
- Measure A = P_max / P_GR at each (σ_w, σ_B); also at multiple t_end values (5, 10, 20).

Parameter range justification (relative to box L=100):

- **σ_w ∈ {2, 5, 10}**: wavepacket width spans 2–10% of the domain. σ_w=2 is the tightest packet for which spectral content from the FFT periodic grid (Δk = 2π/L = 0.063) is still well-sampled (~30 active modes). σ_w=10 is broad enough to overlap a sizeable fraction of σ_B=25 but narrow enough that wavepacket-edge wrap-around is negligible at t_end ≤ 20.
- **σ_B ∈ {10, 25, 50}**: B-field width covers 10–50% of L. σ_B=10 is the marginal case where wavepacket and B-field overlap region is comparable to σ_w; σ_B=25 is the design default (full traversal time ~ 50 units, well-clear of t_end ≤ 20); σ_B=50 approaches the uniform-B₀ limit (sanity check that A converges to Phase B values).
- **t_end ∈ {5, 10, 20}**: brackets the wavepacket transit time σ_w + σ_B. At σ_w + σ_B = 5+10 = 15 the wavepacket has not yet fully cleared the B-field at t_end=10; at σ_w + σ_B = 10 + 50 = 60 the geometry effectively reverts to plane-wave overlap throughout t ≤ 20.

**Quantitative acceptance criterion**: localised geometry is validated for a given (σ_w, σ_B) cell when

```text
A(t_end = 20) / A(t_end = 10) < 1.05
```

i.e. amplification has plateaued to within 5% by t_end=10 (relaxed from the asymptotic limit because we only need t-independent ranking for the inference, not exact convergence). Cells failing this criterion are unusable: the geometry hasn't bounded the growth and we'd be back in the plane-wave-like regime.

Output: `examples/data/v3_localised_geometry_tuning/` with the 9-point × 3-t_end CSV and a 2D heatmap.

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

## Interaction with Phase A-γ (γ_conversion)

Phase A-γ refactors the coupling-aware γ_conversion probe (`log(|target_amp|)/t_test` with multi-`t_test` sampling and log-zero clamp) for a ~1000× amp-chain speedup, deferred to its own session (GH issues #350, #351). Phase E does **not** depend on or block on Phase A-γ:

- Phase E publishes under `P_max:maximize` regardless of A-γ outcome. The geometry pivot is a methodology change (bounded interaction time), independent of which metric drives the likelihood.
- If A-γ ships and validates (Spearman ρ > 0.7 vs P_max), Phase B amp chains *and* Phase E amp chains could both use `gamma_conversion:maximize`. Until then, both use `P_max:maximize`.
- A-γ's `t_test` parameter is unrelated to Phase E's t_end — γ_conversion measures the eigenvalue-implied amplification rate at a probe-internal test time, while Phase E bounds the simulation's interaction time. They're orthogonal.

Decision rule: if A-γ converges *before* Phase E launches, use γ_conversion for Phase E amp chains for the speedup. If not, P_max is the metric. Phase E's geometry conclusions are robust to either choice.

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
- `examples/gertsenshtein/theory_e0_dual_gaussian.toml` — canonical dual-Gaussian theory file (E.0 validated)
- `examples/data/gertsenshtein_e0_dual_gaussian.json` — derived output (6 fields, 22 hamiltonian terms)
- GH #367 — `_evolve_full_matrix` divergence for position-dependent Erf backgrounds
