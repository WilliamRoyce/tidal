# Lagrangian de-pruning audit (v3 architecture, principle 7)

**Created:** 2026-05-10
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)
**Issue:** [#346](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/346)

## Methodology

For each campaign theory, cross-reference:

1. The full quadratic-invariant catalogue in the derivation TOML's `lagrangian` block.
2. The `equations[]` block in the auto-generated JSON.
3. Coefficients fixed to 0 in the campaign script's `--param NAME=0` flags.

Flag any operator that is (a) present in the analytical Lagrangian, (b) **fixed to 0 by analytical-inertness reasoning** (boundary term, total derivative, parity-forbidden under the chosen gauge), and (c) currently excluded from the chain. Per the supervisor's principle 7: re-add such operators with free coefficients; let the chain rediscover their inertness via D_KL ≈ 0 marginal posterior.

Operators fixed to 0 for *non-analytical* reasons (sub-theory definition matching a paper, e.g. "Bahamonde keeps β,δ only") stay fixed — those are conventional, not analytical.

## Audit results

### D1 — Ricci-EM (`torsion_gertsenshtein_nonminimal.json`)

**Free in v2:** `α₁, α₂, α₃, δ₁` (4 params).
**Fixed:** `kappa=1, B0=0.01` (physical constants, not couplings).
**De-pruning candidates:** none. The Ricci-EM theory's Lagrangian has exactly four coupling parameters and all are free.
**v3 action:** rewrite priors (compactified), no parameter set change.

### Stage A — Dark-Photon-Plasma (`dark_photon_plasma.json`)

**Free in v2:** `mA2, deltam, xi, alpha3` (4 params).
**Fixed:** `kappa=1, B0=0.01`.
**De-pruning candidates:** none.
**v3 action:** rewrite priors (compactified), no parameter set change.

### D2.0 Bahamonde — `torsion_gertsenshtein_general_nonminimal.json`

**Free in v2:** `β₁, β₂, β₃, δ₁` (4 params).
**Fixed:** `xi=0, chi=0, zeta1=0, zeta2=0, zeta3=0`.
**De-pruning candidates:** the JSON's Lagrangian includes ξ (kinetic), χ, and ζ_i — fixing them to 0 is a *sub-theory definition* matching Bahamonde et al.'s paper, not an analytical-inertness pruning. Kinetic ξ in particular is necessary to give the torsion-trace propagating dynamics; its absence in Bahamonde is a deliberate choice to study the propagating-torsion-mass-only regime.
**v3 action:** keep the sub-theory partition (it's a literature-comparison choice). Rewrite priors only.

### D2.1 Barker — same JSON as D2.0

**Free in v2:** `β₁, β₂, β₃, ξ, χ` (5 params).
**Fixed:** `delta1=0, zeta1=0, zeta2=0, zeta3=0`.
**De-pruning candidates:** Barker's paper omits δ₁ and ζ_i to focus on the χ-only nonminimal coupling. Same sub-theory-definition rationale as D2.0.
**v3 action:** keep partition; rewrite priors only.

### D2.2 Shapiro — same JSON

**Free in v2:** `β₁, β₂, β₃, ζ₁, ζ₂, ζ₃` (6 params).
**Fixed:** `xi=0, delta1=0, chi=0`.
**De-pruning candidates:** Shapiro's paper studies derivative-torsion couplings ζ_i without ξ/δ₁/χ. Sub-theory.
**v3 action:** keep partition; rewrite priors only.

### D2.3 Full T5 — same JSON, 9 params free

**Free in v2:** `β₁, β₂, β₃, ξ, δ₁, χ, ζ₁, ζ₂, ζ₃` (9 params; the full Lagrangian).
**Fixed:** none beyond physical constants.
**De-pruning candidates:** none (already maximal).
**v3 action:** rewrite priors only.

### T6 parity-odd YM-PGT — `torsion_gertsenshtein_parity_odd.json` (NOT in Phase A; deferred)

**Free (when run):** T5 (9) + d14, d15, d16, d17, d19, d20, d21, zt1–zt6 = 22 params per the supervisor's May 8 note.
**Lagrangian:** the JSON's `lagrangian_expr` already includes all the parity-odd torsion-quadratic (d14–17), parity-odd Ricci/Riemann–F (d19–21), and parity-odd ∇T–F (zt1–6) terms.
**Possible de-pruning candidate:** the parity-odd FF̃ term `ε^{abcd} F_{ab} F_{cd}` (the Pontryagin density of EM, a topological / total-derivative invariant). This term is **not currently in the JSON's Lagrangian** because it analytically integrates to a boundary and contributes nothing to the EOM in flat spacetime. By the supervisor's principle 7, T6 should include it with a free coefficient (call it `d18` to fill the gap in the d14–21 enumeration) so the chain can rediscover its inertness via D_KL ≈ 0.
**v3 action:** **deferred** — Phase A's 6 campaigns don't include T6. When T6 is brought into the v3 campaign queue (post-Phase B), re-derive the JSON with the FF̃ term added; the v3 likelihood architecture handles the wider prior gracefully.

### EH — Einstein-Maxwell + Euler-Heisenberg (NOT in Phase A; deferred)

**Status:** derivation pending (Wolfram-side issue per supervisor meeting note §"Roadmap"). Once derived, audit will check whether any F⁴ contractions were dropped on analytical grounds (some `(F·F)·(F̃·F̃)` terms can vanish identically). Action deferred until the JSON exists.

## Summary

For Phase A's 6 campaigns (D1, Stage A, D2.0–D2.3), **no analytical-inertness de-pruning is needed** — all currently-pruned parameters are pruned for sub-theory-definition or sampling-strategy reasons, not because they're known boundary/topological invariants. The principle-7 application is concentrated in T6 (parity-odd, deferred to a post-Phase-B campaign) and EH (deferred until derived).

This audit is therefore mostly a *negative finding* for Phase A: the v3_permissive campaign scripts can rewrite priors without touching the parameter-set partitions established by the v2 sub-theory definitions. The principle is logged and will activate when T6 / EH come online.

## References

- [docs/V3_ARCHITECTURE.md](V3_ARCHITECTURE.md) §"Why v3?" #7
- [docs/meetings/2026-05-08_supervisor.md](meetings/2026-05-08_supervisor.md)
- Original D2 sub-theory script: [scripts/campaign/submit_campaign_D.sh](../scripts/campaign/submit_campaign_D.sh)
- Bahamonde et al. (2024); Barker (arXiv:2406.12826); Shapiro (arXiv:hep-th/0103093) — define the sub-theory partitions.
