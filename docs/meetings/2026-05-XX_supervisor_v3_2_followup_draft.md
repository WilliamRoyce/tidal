# Supervisor meeting draft — v3.2 prior/ghost follow-up

**Status:** DRAFT — to be sent in next regular meeting / email
**Created:** 2026-05-11
**Companion to:** [docs/V3_2_DESIGN_INVESTIGATION.md](../V3_2_DESIGN_INVESTIGATION.md), [docs/V3_ARCHITECTURE.md](../V3_ARCHITECTURE.md)

## Subject

v3.2: kinetic vs mass-like parameter classification (and PSALTer integration timeline)

## What we've done since 8 May

Phase B v3.1 chains are landing. Headline so far (all under v3 architecture: tachyon-permissive, soft floor, compactified arctan_uniform on sign-symmetric dims, log_uniform on positive-only kinetic):

| Chain | v2 reference | v2 log Z | v3 log Z | v3 joint D_KL |
|---|---|---|---|---|
| D1 amp smoke (29149987) | 28896653 | +0.72 | **+13.29** | **4.77 nats** |
| Stage A amp (29189966) | 28474676 | −0.07 | **+9.31** | **6.59 nats** |
| Stage A sup (29199129) | 28477675 | +0.65 | **+4.08** | **7.30 nats** |

Per-coupling MAPs in v3.1 land outside v2's narrow priors in several places (α₁/δ₁ in D1, mA²/deltam/α₃/xi in Stage A). v2's "Stage A null verdict" turns out to have been a prior-support artefact: the chain wanted higher mA² and larger deltam than v2 allowed, and v3's wider priors find clear signal there. Detail and corner plots in `docs/comparison/`.

## What we're stuck on

A v2-inherited assumption that the v3 architecture should challenge but currently doesn't.

**Concretely**: we still treat α₃ (D1), and mA² / ξ / α₃ (Stage A), as positive-only via `log_uniform:1e-3:1e3`. But:

- v3 is tachyon-permissive — we *allow* m² < 0 (the dark-photon tachyonic regime is a legitimate parameter region under v3, not a forbidden one).
- For mA² in particular, which appears in the Stage A Lagrangian as a mass-squared term `mA²·A_μA^μ`, sampling it as positive-only assumes "no tachyonic dark photon" — exactly what v3 is supposed to let the chain decide.
- The Stage A sup chain (29199129) just landed with α₃ MAP = 0.00117, hitting our prior lower bound `1e-3`. The chain wants to push α₃ smaller (or possibly negative, if α₃ is mass-like in this theory).

The thing that should be positive-only by physics is the **kinetic coefficient** (negative kinetic → ghost → Hamiltonian unbounded below → quantum-mechanical pathology). The thing that should be sign-symmetric in v3 is everything **mass-like**.

## Our specific question

For each parameter currently sampled as `log_uniform:1e-3:1e3` in v3.1, is it:

(a) **Kinetic-coefficient** — multiplies a `(∂φ)²` term, positivity required to avoid ghosts — or

(b) **Mass-like** — multiplies a `φ²` algebraic term, sign-symmetric under v3 tachyon-permissive policy?

Specifically:

| Theory | Parameter | What is it physically? |
|---|---|---|
| D1 (torsion_gertsenshtein_nonminimal) | α₃ | One of the dimensionless torsion-quadratic invariant coefficients — kinetic or algebraic? |
| Stage A (dark_photon_plasma) | mA² | Dark-photon mass-squared — confirm: this should be sign-symmetric under v3, right? |
| Stage A | ξ | Kinetic / mass-mixing / dimensionless coupling? |
| Stage A | α₃ | Torsion-photon coupling — same operator family as D1's α₃? |
| D2.0–D2.3 sub-theories | (various) | We'll inventory and follow up if any are currently `log_uniform` |

Once classified, we can switch the mass-like ones to `arctan_uniform:-89:89` and rerun the affected chains. The kinetic-only ones stay `log_uniform` until v3.2 (per-coupling abs-fold in the cubed-sphere joint prior) lands.

## Why we're not just doing it ourselves

Your point on PSALTer integration in the 8 May meeting still applies: a homegrown ghost check (just checking the Hamiltonian density for negative kinetic terms at each chain sample) is unreliable for PGT theories because gauge-fixing and constraint elimination mix field components — saturated-propagator residue analysis is what's actually needed, and PSALTer v2's kinetic-matrix definiteness approach is the correct algorithm.

We've filed an issue tracking PSALTer integration as the ultimate goal (Phase C+ timeline, gated on PSALTer's Python loader being public API). Until then, your classification answers above let us apply the lighter-touch "per-coupling abs-fold" interim — we get the joint prior covering all couplings uniformly without committing to the full PSALTer bridge yet.

## Cross-references for you

- `docs/V3_ARCHITECTURE.md` — canonical v3 architecture, soft-penalty table, compactified priors
- `docs/V3_2_DESIGN_INVESTIGATION.md` — full investigation notes
- `docs/comparison/*` — per-chain v2 vs v3 numbers and corner plots
- GH `[v3.2-classify]` issue — this question, tracked
- GH `[v3.2-prior]` issue — the per-coupling abs-fold implementation, gated on this
- GH `[v3.2-psalter]` issue — the long-term integration target
