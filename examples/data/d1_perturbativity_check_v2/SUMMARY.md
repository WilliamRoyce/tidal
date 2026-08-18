# Phase 6.L — Perturbativity validation at v2 amp MAPs

**Date:** 2026-05-08
**Branch:** `hpc/pgt-survey`
**Probe HEAD:** `98c87d7` (τ=0.15)

## Summary

All three v2 amp MAPs (τ=0.15 probe) sit at the **boundary attractor**
γ_eff ≈ 0.13–0.15. The probe correctly admits them, but the t_end
sweep reveals **slow tachyonic growth of A within the probe pass band**:
A(20)/A(10) ≈ 2.4–3.2 across all three MAPs (FAIL the [0.5, 2.0]
perturbativity band).

## Per-MAP results (B₀ = 1e-4 throughout)

| MAP | Job | Grid/L | γ_eff (probe) | A(t=10) | A(t=25) | A(20)/A(10) | log A slope |
|-----|-----|--------|---------------|---------|---------|-------------|-------------|
| 28896653 | hi-res | 512/100 | **0.1367** | 38.85 | 254.08 | **3.17** | 0.1060/unit |
| 28883112 | grid=256 | 256/100 | **0.1500** (boundary) | 20.00 | 85.76 | **2.40** | 0.0847/unit |
| 28982029 | INTR-reduced | 256/100 | **0.1500** (boundary) | (similar) | (similar) | **~2.4** | ~0.085/unit |

Raw P growth (corrected for sin² baseline scaling): **log P slope = 0.24–0.26 per unit** at all three MAPs — i.e., P ∝ exp(0.25·t) approximately, indicating effective field growth γ_h ≈ 0.12–0.13.

## Interpretation

The **τ=0.15 probe correctly bounds the worst-case growth at exp(0.15·t_test)
= exp(3) ≈ 20× over t_test=20**, but it does NOT guarantee the [0.5, 2.0]
perturbativity pass band over t_end=20 vs t_end=10. The chain finds
the maximum-A point under the probe constraint, which by construction
sits at the boundary where γ ≈ 0.13–0.15.

**Publication implications:** the A_max ≈ 38 reported at t_end=10 by
the chain is **partially a slow tachyonic enhancement** that does not
survive to t_end=25. The t-independent "static" component (light-mediator
enhancement à la 1/m²) is approximately A_static ≈ A(t=10) / exp(γ·t_end)
≈ 38 / exp(1.37) ≈ **10**.

## Probe-truncation artifact

An earlier check at the **3-decimal truncated** MAP coords gave γ_eff =
0.0001 for 28896653 — orders of magnitude below the full-precision
result. This is because the chain MAP sits in a **very narrow ridge**
of stability (the boundary attractor); rounding the parameters by
0.001 moves off the ridge into a different stability regime. Always
use **full-precision MAP coordinates** for probe sanity checks.

## Action items for #340 / #341

The τ=0.15 probe does its job (caps growth at exp(0.15·t_test) ≈ 20×),
but the chain MAP attractor sits exactly at this boundary, so the
"publication A_max" includes residual t_end-dependent growth.

**Two options for the manuscript:**

1. **Tighten further to τ=0.07** (or similar). Caps growth at
   exp(0.07·10) ≈ 2× — properly perturbative. Will shift the chain
   log Z (likely lower) and steer the MAP to a less-extremal point.
2. **Report A(t_end) at multiple t_end values** for transparency.
   Disclose that A_max=38 at t_end=10 includes ~4× slow-growth
   factor; A_static ≈ 10.

Recommended: try τ=0.07 in a single test INTR run and see what log Z
results; if not too costly, adopt as the canonical probe.

## Files

- `28896653_b0_sweep/results.csv` — B₀ sweep (PASS, 1.5% variation)
- `28896653_tend_sweep/t{5,10,15,20,25}/` — t_end sweep (FAIL)
- `28883112_tend_sweep/t{5,10,15,20,25}/` — t_end sweep (FAIL)
